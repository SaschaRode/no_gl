#include "no_gl.h"

#define MAX_JOBS 128
#define MAX_TEXTURES 8
#define MAX_UNIFORMS 32
#define MAX_TRIANGLES (32*4096)

typedef void *thread_callback_t(void *data);

typedef struct thread_t
{
    uint64_t opaque;
} thread_t;

typedef struct semaphore_t
{
    uint64_t opaque[4];
} semaphore_t;

typedef void job_callback_t(void *data);

typedef struct job_t
{
    job_callback_t *callback;
    void *data;
} job_t;

typedef struct job_queue_t
{
    uint32_t num_threads;
    int32_t volatile active_jobs;
    int32_t volatile current_write;
    int32_t volatile current_read;
    job_t jobs[MAX_JOBS];
    semaphore_t sem;
} job_queue_t;

typedef struct rect_t
{
    int32_t x, y, w, h;
} rect_t;

typedef union uniform_t
{
    scalar_t scalar;
    svec2_t vec2;
    svec3_t vec3;
    svec4_t vec4;
} uniform_t;

typedef struct homogeneous_vertex_t
{
    vec4_t hpos;
    vec3_t pos;
    vec3_t normal;
    vec3_t tangent;
    vec2_t uv;
} homogeneous_vertex_t;

typedef struct triangle_t
{
    homogeneous_vertex_t v0, v1, v2;
} triangle_t;

typedef struct no_gl_t
{
    framebuffer_t *framebuffer;
    fragment_shader_t *fragment_shader;
    texture_t *texture[MAX_TEXTURES];
    uniform_t uniform[MAX_UNIFORMS];
    triangle_t triangles[MAX_TRIANGLES];
    uint32_t num_triangles;
    job_queue_t job_queue;
    m4x4_t model;
    m4x4_t view;
    m4x4_t proj;
    xm4x4_t xmodel;
    xm4x4_t xview;
    xm4x4_t xproj;
    bool valid;
} no_gl_t;

static no_gl_t gl = { 0 };

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#define READ_BARRIER _ReadBarrier()
#define WRITE_BARRIER _WriteBarrier()

static int32_t atomic_add32(int32_t volatile *value, int32_t addend)
{
    int32_t result = _InterlockedExchangeAdd((long volatile *)value, addend);
    return(result);
}

static int32_t atomic_compare_swap32(int32_t volatile *value, int32_t exch, int32_t comp)
{
    int32_t result = _InterlockedCompareExchange((long volatile *)value, exch, comp);
    return(result);
}

static semaphore_t create_semaphore(uint32_t value)
{
    semaphore_t sem =
        {
            (uint64_t)CreateSemaphoreEx(0, value, INT16_MAX, 0, 0, SEMAPHORE_ALL_ACCESS),
        };
    return(sem);
}

static void post_semaphore(semaphore_t *sem)
{
    ReleaseSemaphore((HANDLE)sem->opaque[0], 1, 0);
}

static void wait_for_semaphore(semaphore_t *sem)
{
    WaitForSingleObjectEx((HANDLE)sem->opaque[0], INFINITE, false);
}

static void destroy_semaphore(semaphore_t *sem)
{
    CloseHandle((HANDLE)sem->opaque[0]);
}

static thread_t create_thread(thread_callback_t *callback, void *data)
{
    thread_t thread;
    typedef DWORD win32_thread_callback_t(LPVOID Parameter);
    thread.opaque = (uint64_t)CreateThread(0, 0, (win32_thread_callback_t *)callback, data, 0, 0);
    return(thread);
}

static void destroy_thread(thread_t thread)
{
    CloseHandle((HANDLE)thread.opaque);
}

static void wait_for_thread(thread_t thread)
{
    WaitForSingleObject((HANDLE)thread.opaque, INFINITE);
}

static uint32_t get_num_processors(void)
{
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return(info.dwNumberOfProcessors);
}
#endif

static void *job_queue_thread(void *data)
{
    job_queue_t *queue = data;

    for(;;)
    {
        int32_t read = queue->current_read;
        int32_t next_read = (read + 1)%MAX_JOBS;

        if(read != queue->current_write)
        {
            int32_t result = atomic_compare_swap32(&queue->current_read, next_read, read);
            if(result == read)
            {
                job_t *job = queue->jobs + read;
                job->callback(job->data);
                atomic_add32(&queue->active_jobs, -1);
            }
        }
        else
        {
            wait_for_semaphore(&queue->sem);
        }
    }
}

static void wait_for_all_jobs(void)
{
    while(gl.job_queue.active_jobs != 0)
    {
    }
}

static void add_job(job_callback_t *callback, void *data)
{
    job_t *job = gl.job_queue.jobs + gl.job_queue.current_write;
    job->callback = callback;
    job->data = data;
    uint32_t next_write = (gl.job_queue.current_write + 1)%MAX_JOBS;
    WRITE_BARRIER;
    atomic_add32(&gl.job_queue.active_jobs, 1);
    gl.job_queue.current_write = next_write;
    post_semaphore(&gl.job_queue.sem);
}

static void setup_job_queue(void)
{
    gl.job_queue.num_threads = get_num_processors();
    gl.job_queue.sem = create_semaphore(0);
    for(uint32_t i = 0; i < gl.job_queue.num_threads; i++)
        create_thread(&job_queue_thread, &gl.job_queue);
}

static void set_depth(framebuffer_t *framebuffer, int32_t x, int32_t y, float depth)
{
    framebuffer->depth[y*framebuffer->w + x] = depth;
}

static void set_color(framebuffer_t *framebuffer, int32_t x, int32_t y, uint32_t color)
{
    framebuffer->color[y*framebuffer->w + x] = color;
}

static float get_depth(framebuffer_t *framebuffer, int32_t x, int32_t y)
{
    float depth = framebuffer->depth[y*framebuffer->w + x];
    return(depth);
}

static uint32_t get_color(framebuffer_t *framebuffer, int32_t x, int32_t y)
{
    uint32_t color = framebuffer->color[y*framebuffer->w + x];
    return(color);
}

static uint32_t pack_color(vec4_t color)
{
    uint32_t color32 = 0;
    color32 |= (uint8_t)(clamp01(color.b)*255.0f);
    color32 |= (uint8_t)(clamp01(color.g)*255.0f) << 8;
    color32 |= (uint8_t)(clamp01(color.r)*255.0f) << 16;
    color32 |= (uint8_t)(clamp01(color.a)*255.0f) << 24;
    return(color32);
}

static vec4_t unpack_color(uint32_t color32)
{
    float inv_255 = 1.0f/255.0f;
    vec4_t color =
        {
            (float)((color32 >> 16) & 0xFF)*inv_255,
            (float)((color32 >> 8) & 0xFF)*inv_255,
            (float)(color32 & 0xFF)*inv_255,
            (float)((color32 >> 24) & 0xFF)*inv_255,
        };
    return(color);
}

static __m128i xpack_color(xvec4_t color)
{
    __m128 min_value = _mm_set_ps1(0.0f);
    __m128 max_value = _mm_set_ps1(255.0f);

    __m128 b = _mm_max_ps(_mm_min_ps(_mm_mul_ps(color.b, max_value), max_value), min_value);
    __m128i color32 = _mm_cvttps_epi32(b);
    __m128 g = _mm_max_ps(_mm_min_ps(_mm_mul_ps(color.g, max_value), max_value), min_value);
    color32 = _mm_or_si128(color32, _mm_slli_epi32(_mm_cvttps_epi32(g), 8));
    __m128 r = _mm_max_ps(_mm_min_ps(_mm_mul_ps(color.r, max_value), max_value), min_value);
    color32 = _mm_or_si128(color32, _mm_slli_epi32(_mm_cvttps_epi32(r), 16));
    __m128 a = _mm_max_ps(_mm_min_ps(_mm_mul_ps(color.a, max_value), max_value), min_value);
    color32 = _mm_or_si128(color32, _mm_slli_epi32(_mm_cvttps_epi32(a), 24));

    return(color32);
}

static xvec4_t xunpack_color(__m128i color32)
{
    __m128i color_mask = _mm_set1_epi32(0xFF);
    __m128 inv_255 = _mm_set1_ps(1.0f/255.0f);

    __m128i r = _mm_and_si128(_mm_srli_epi32(color32, 16), color_mask);
    __m128i g = _mm_and_si128(_mm_srli_epi32(color32, 8), color_mask);
    __m128i b = _mm_and_si128(color32, color_mask);
    __m128i a = _mm_and_si128(_mm_srli_epi32(color32, 24), color_mask);

    xvec4_t color =
        {
            _mm_mul_ps(_mm_cvtepi32_ps(r), inv_255),
            _mm_mul_ps(_mm_cvtepi32_ps(g), inv_255),
            _mm_mul_ps(_mm_cvtepi32_ps(b), inv_255),
            _mm_mul_ps(_mm_cvtepi32_ps(a), inv_255),
        };
    return(color);
}

static homogeneous_vertex_t clip_edge(homogeneous_vertex_t *v0, homogeneous_vertex_t *v1)
{
    float t = (-v1->hpos.w - v1->hpos.z) / (v0->hpos.z + v0->hpos.w - v1->hpos.w - v1->hpos.z);
    homogeneous_vertex_t result =
        {
            .hpos = vec4_lerp(v1->hpos, v0->hpos, t),
            .pos = vec3_lerp(v1->pos, v0->pos, t),
            .normal = vec3_lerp(v1->normal, v0->normal, t),
            .tangent = vec3_lerp(v1->tangent, v0->tangent, t),
            .uv = vec2_lerp(v1->uv, v0->uv, t),
        };
    return(result);
}

static void convert_vertex(vertex_t *v, homogeneous_vertex_t *result)
{
    result->pos = v->pos;
    result->normal = v->normal;
    result->tangent = v->tangent;
    result->uv = v->uv;
}

static float edge_func(vec2_t a, vec2_t b, vec2_t c)
{
    return((c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x));
}

static __m128 xedge_func(__m128 ax, __m128 ay, __m128 bx, __m128 by, __m128 cx, __m128 cy)
{
    return(_mm_sub_ps(_mm_mul_ps(_mm_sub_ps(cx, ax), _mm_sub_ps(by, ay)), _mm_mul_ps(_mm_sub_ps(cy, ay), _mm_sub_ps(bx, ax))));
}

static void draw_region(rect_t *region)
{
#if defined(FAST_MODE)
    __m128 zero = _mm_set_ps1(0.0f);
    __m128 one = _mm_set_ps1(1.0f);

    __m128 offset_x = _mm_set_ps(3.5f, 2.5f, 1.5f, 0.5f);
    __m128 offset_y = _mm_set_ps1(0.5f);

    for(uint32_t i = 0; i < gl.num_triangles; i++)
    {
        triangle_t *triangle = &gl.triangles[i];

        homogeneous_vertex_t *v0 = &triangle->v0;
        homogeneous_vertex_t *v1 = &triangle->v1;
        homogeneous_vertex_t *v2 = &triangle->v2;

        vec2_t bounds_min = vec2_min(vec2_min(v0->hpos.xy, v1->hpos.xy), v2->hpos.xy);
        vec2_t bounds_max = vec2_max(vec2_max(v0->hpos.xy, v1->hpos.xy), v2->hpos.xy);

        // TODO: handle screen sizes that are not a multiple of four.
        int32_t start_x = MAX((int32_t)(floorf(bounds_min.x/4.0f)*4.0f), region->x);
        int32_t start_y = MAX((int32_t)floorf(bounds_min.y), region->y);
        int32_t end_x = MIN((int32_t)(ceilf(bounds_max.x/4.0f)*4.0f), (region->x + region->w));
        int32_t end_y = MIN((int32_t)ceilf(bounds_max.y), (region->y + region->h));

        if(start_x >= end_x || start_y >= end_y)
            continue;

        xvec3_t p0 = xvec3_setv(vec3_muls(v0->pos, v0->hpos.w));
        xvec3_t p1 = xvec3_setv(vec3_muls(v1->pos, v1->hpos.w));
        xvec3_t p2 = xvec3_setv(vec3_muls(v2->pos, v2->hpos.w));

        xvec3_t n0 = xvec3_setv(vec3_muls(v0->normal, v0->hpos.w));
        xvec3_t n1 = xvec3_setv(vec3_muls(v1->normal, v1->hpos.w));
        xvec3_t n2 = xvec3_setv(vec3_muls(v2->normal, v2->hpos.w));

        xvec3_t t0 = xvec3_setv(vec3_muls(v0->tangent, v0->hpos.w));
        xvec3_t t1 = xvec3_setv(vec3_muls(v1->tangent, v1->hpos.w));
        xvec3_t t2 = xvec3_setv(vec3_muls(v2->tangent, v2->hpos.w));

        xvec2_t uv0 = xvec2_setv(vec2_muls(v0->uv, v0->hpos.w));
        xvec2_t uv1 = xvec2_setv(vec2_muls(v1->uv, v1->hpos.w));
        xvec2_t uv2 = xvec2_setv(vec2_muls(v2->uv, v2->hpos.w));

        xvec4_t hp0 = xvec4_setv(v0->hpos);
        xvec4_t hp1 = xvec4_setv(v1->hpos);
        xvec4_t hp2 = xvec4_setv(v2->hpos);

        float area = edge_func(v0->hpos.xy, v1->hpos.xy, v2->hpos.xy);
        __m128 inv_area = _mm_set_ps1(1.0f/area);

        for(int32_t y = start_y; y < end_y; y++)
        {
            for(int32_t x = start_x; x < end_x; x += 4)
            {
                uint32_t *color_ptr = &gl.framebuffer->color[y*gl.framebuffer->w + x];
                _mm_prefetch((char *)color_ptr, _MM_HINT_T0);

                float *depth_ptr = &gl.framebuffer->depth[y*gl.framebuffer->w + x];
                _mm_prefetch((char *)depth_ptr, _MM_HINT_T0);

                __m128 px = _mm_add_ps(_mm_set_ps1((float)x), offset_x);
                __m128 py = _mm_add_ps(_mm_set_ps1((float)y), offset_y);

                __m128 w0 = xedge_func(hp1.x, hp1.y, hp2.x, hp2.y, px, py);
                __m128 w1 = xedge_func(hp2.x, hp2.y, hp0.x, hp0.y, px, py);
                __m128 w2 = xedge_func(hp0.x, hp0.y, hp1.x, hp1.y, px, py);

                __m128 inside_compare = _mm_and_ps(_mm_and_ps(_mm_cmple_ps(w0, zero), _mm_cmple_ps(w1, zero)), _mm_cmple_ps(w2, zero));

                if(_mm_testz_si128(_mm_castps_si128(inside_compare), _mm_castps_si128(inside_compare)) != 0)
                    continue;

                w0 = _mm_mul_ps(w0, inv_area);
                w1 = _mm_mul_ps(w1, inv_area);
                w2 = _mm_mul_ps(w2, inv_area);

                __m128 dst_depth = _mm_add_ps(_mm_add_ps(_mm_mul_ps(w0, hp0.z), _mm_mul_ps(w1, hp1.z)), _mm_mul_ps(w2, hp2.z));
                __m128 src_depth = _mm_load_ps(depth_ptr);
                __m128 depth_compare = _mm_cmplt_ps(dst_depth, src_depth);

                if(_mm_testz_si128(_mm_castps_si128(depth_compare), _mm_castps_si128(depth_compare)) != 0)
                    continue;

                __m128 w = _mm_add_ps(_mm_add_ps(_mm_mul_ps(w0, hp0.w), _mm_mul_ps(w1, hp1.w)), _mm_mul_ps(w2, hp2.w));
                __m128 inv_w = _mm_div_ps(one, w);

                xvec3_t pos;
                pos.x = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(p0.x, w0), _mm_mul_ps(p1.x, w1)), _mm_mul_ps(p2.x, w2)), inv_w);
                pos.y = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(p0.y, w0), _mm_mul_ps(p1.y, w1)), _mm_mul_ps(p2.y, w2)), inv_w);
                pos.z = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(p0.z, w0), _mm_mul_ps(p1.z, w1)), _mm_mul_ps(p2.z, w2)), inv_w);

                xvec3_t normal;
                normal.x = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(n0.x, w0), _mm_mul_ps(n1.x, w1)), _mm_mul_ps(n2.x, w2)), inv_w);
                normal.y = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(n0.y, w0), _mm_mul_ps(n1.y, w1)), _mm_mul_ps(n2.y, w2)), inv_w);
                normal.z = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(n0.z, w0), _mm_mul_ps(n1.z, w1)), _mm_mul_ps(n2.z, w2)), inv_w);

                xvec3_t tangent;
                tangent.x = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(t0.x, w0), _mm_mul_ps(t1.x, w1)), _mm_mul_ps(t2.x, w2)), inv_w);
                tangent.y = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(t0.y, w0), _mm_mul_ps(t1.y, w1)), _mm_mul_ps(t2.y, w2)), inv_w);
                tangent.z = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(t0.z, w0), _mm_mul_ps(t1.z, w1)), _mm_mul_ps(t2.z, w2)), inv_w);

                xvec2_t uv;
                uv.x = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(uv0.x, w0), _mm_mul_ps(uv1.x, w1)), _mm_mul_ps(uv2.x, w2)), inv_w);
                uv.y = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(uv0.y, w0), _mm_mul_ps(uv1.y, w1)), _mm_mul_ps(uv2.y, w2)), inv_w);

                xvec4_t dst_color = gl.fragment_shader(pos, normal, tangent, uv);
                xvec4_t src_color = xunpack_color(_mm_load_si128((__m128i *)color_ptr));
                dst_color.rgb = xvec3_lerp(src_color.rgb, dst_color.rgb, dst_color.a);

                __m128i mask = _mm_castps_si128(_mm_and_ps(inside_compare, depth_compare));
                _mm_maskstore_epi32((int32_t *)color_ptr, mask, xpack_color(dst_color));
                _mm_maskstore_ps(depth_ptr, mask, dst_depth);
            }
        }
    }
#else
    for(uint32_t i = 0; i < gl.num_triangles; i++)
    {
        triangle_t *triangle = &gl.triangles[i];

        homogeneous_vertex_t *v0 = &triangle->v0;
        homogeneous_vertex_t *v1 = &triangle->v1;
        homogeneous_vertex_t *v2 = &triangle->v2;

        vec2_t bounds_min = vec2_min(vec2_min(v0->hpos.xy, v1->hpos.xy), v2->hpos.xy);
        vec2_t bounds_max = vec2_max(vec2_max(v0->hpos.xy, v1->hpos.xy), v2->hpos.xy);

        int32_t start_x = MAX((int32_t)floorf(bounds_min.x), region->x);
        int32_t start_y = MAX((int32_t)floorf(bounds_min.y), region->y);
        int32_t end_x = MIN((int32_t)ceilf(bounds_max.x), (region->x + region->w));
        int32_t end_y = MIN((int32_t)ceilf(bounds_max.y), (region->y + region->h));

        if(start_x >= end_x || start_y >= end_y)
            continue;

        vec3_t p0 = vec3_muls(v0->pos, v0->hpos.w);
        vec3_t p1 = vec3_muls(v1->pos, v1->hpos.w);
        vec3_t p2 = vec3_muls(v2->pos, v2->hpos.w);

        vec3_t n0 = vec3_muls(v0->normal, v0->hpos.w);
        vec3_t n1 = vec3_muls(v1->normal, v1->hpos.w);
        vec3_t n2 = vec3_muls(v2->normal, v2->hpos.w);

        vec3_t t0 = vec3_muls(v0->tangent, v0->hpos.w);
        vec3_t t1 = vec3_muls(v1->tangent, v1->hpos.w);
        vec3_t t2 = vec3_muls(v2->tangent, v2->hpos.w);

        vec2_t uv0 = vec2_muls(v0->uv, v0->hpos.w);
        vec2_t uv1 = vec2_muls(v1->uv, v1->hpos.w);
        vec2_t uv2 = vec2_muls(v2->uv, v2->hpos.w);

        float area = edge_func(v0->hpos.xy, v1->hpos.xy, v2->hpos.xy);
        float inv_area = 1.0f/area;

        for(int32_t y = start_y; y < end_y; y++)
        {
            for(int32_t x = start_x; x < end_x; x++)
            {
                vec2_t pixel = { (float)x + 0.5f, (float)y + 0.5f };

                float w0 = edge_func(v1->hpos.xy, v2->hpos.xy, pixel);
                float w1 = edge_func(v2->hpos.xy, v0->hpos.xy, pixel);
                float w2 = edge_func(v0->hpos.xy, v1->hpos.xy, pixel);

                if(w0 <= 0 && w1 <= 0 && w2 <= 0)
                {
                    w0 *= inv_area;
                    w1 *= inv_area;
                    w2 *= inv_area;

                    float depth = (w0*v0->hpos.z + w1*v1->hpos.z + w2*v2->hpos.z);
                    if(depth < get_depth(gl.framebuffer, x, y))
                    {
                        float inv_w = 1.0f/(w0*v0->hpos.w + w1*v1->hpos.w + w2*v2->hpos.w);

                        vec3_t pos = vec3_muls(vec3_add(vec3_add(vec3_muls(p0, w0), vec3_muls(p1, w1)), vec3_muls(p2, w2)), inv_w);
                        vec3_t normal = vec3_muls(vec3_add(vec3_add(vec3_muls(n0, w0), vec3_muls(n1, w1)), vec3_muls(n2, w2)), inv_w);
                        vec3_t tangent = vec3_muls(vec3_add(vec3_add(vec3_muls(t0, w0), vec3_muls(t1, w1)), vec3_muls(t2, w2)), inv_w);
                        vec2_t uv = vec2_muls(vec2_add(vec2_add(vec2_muls(uv0, w0), vec2_muls(uv1, w1)), vec2_muls(uv2, w2)), inv_w);

                        vec4_t dst_color = gl.fragment_shader(pos, normal, tangent, uv);
                        vec4_t src_color = unpack_color(get_color(gl.framebuffer, x, y));
                        dst_color.rgb = vec3_lerp(src_color.rgb, dst_color.rgb, dst_color.a);

                        set_color(gl.framebuffer, x, y, pack_color(dst_color));
                        set_depth(gl.framebuffer, x, y, depth);
                    }
                }
            }
        }
    }
#endif
}

framebuffer_t *create_framebuffer(int32_t w, int32_t h)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(w > 0 && (w % 4) == 0); // framebuffer width has to be greater zero and multiple of four!
    ASSERT(h > 0 && (h % 4) == 0); // framebuffer height has to be greater zero and multiple of four!

    framebuffer_t *framebuffer = malloc(sizeof(framebuffer_t));
    *framebuffer = (framebuffer_t)
        {
            .w = w,
            .h = h,
            .color = malloc(w*h*sizeof(uint32_t)),
            .depth = malloc(w*h*sizeof(float)),
        };
    return(framebuffer);
}

texture_t *create_texture(int32_t w, int32_t h, address_mode_t mode)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(w > 0 && w == next_pow2(w)); // texture width has to be greater zero and power of two!
    ASSERT(h > 0 && h == next_pow2(h)); // texture height has to be greater zero and power of two!

    texture_t *texture = malloc(sizeof(texture_t));
    *texture = (texture_t)
        {
            .w = w,
            .h = h,
            .pixels = malloc(w*h*sizeof(uint32_t)),
            .address_mode = mode,
        };
    return(texture);
}

void free_framebuffer_memory(framebuffer_t *framebuffer)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(framebuffer);

    free(framebuffer->color);
    free(framebuffer->depth);
    free(framebuffer);
}

void free_texture_memory(texture_t *texture)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(texture);

    free(texture->pixels);
    free(texture);
}

void set_model_matrix(m4x4_t matrix)
{
    ASSERT(gl.valid); // call setup_gl() first!

    gl.model = m4x4_transpose(matrix);
    gl.xmodel = xm4x4_setm(matrix);
}

void set_view_matrix(m4x4_t matrix)
{
    ASSERT(gl.valid); // call setup_gl() first!

    gl.view = matrix;
    gl.xview = xm4x4_setm(matrix);
}

void set_projection_matrix(m4x4_t matrix)
{
    ASSERT(gl.valid); // call setup_gl() first!

    gl.proj = matrix;
    gl.xproj = xm4x4_setm(matrix);
}

void set_texture(uint32_t idx, texture_t *texture)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_TEXTURES);

    gl.texture[idx] = texture;
}

texture_t *get_texture(uint32_t idx)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_TEXTURES);

    return(gl.texture[idx]);
}

sm4x4_t get_model_matrix(void)
{
    ASSERT(gl.valid); // call setup_gl() first!

#if defined(FAST_MODE)
    return(gl.xmodel);
#else
    return(gl.model);
#endif
}

sm4x4_t get_view_matrix(void)
{
    ASSERT(gl.valid); // call setup_gl() first!

#if defined(FAST_MODE)
    return(gl.xview);
#else
    return(gl.view);
#endif
}

sm4x4_t get_projection_matrix(void)
{
    ASSERT(gl.valid); // call setup_gl() first!

#if defined(FAST_MODE)
    return(gl.xproj);
#else
    return(gl.proj);
#endif
}

void set_uniform_scalar(uint32_t idx, float value)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

    gl.uniform[idx].scalar = scalar(value);
}

void set_uniform_vec2(uint32_t idx, vec2_t value)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

#if defined FAST_MODE
    gl.uniform[idx].vec2 = xvec2_setv(value);
#else
    gl.uniform[idx].vec2 = value;
#endif
}

void set_uniform_vec3(uint32_t idx, vec3_t value)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

#if defined FAST_MODE
    gl.uniform[idx].vec3 = xvec3_setv(value);
#else
    gl.uniform[idx].vec3 = value;
#endif
}

void set_uniform_vec4(uint32_t idx, vec4_t value)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

#if defined FAST_MODE
    gl.uniform[idx].vec4 = xvec4_setv(value);
#else
    gl.uniform[idx].vec4 = value;
#endif
}

scalar_t get_uniform_scalar(uint32_t idx)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

    return(gl.uniform[idx].scalar);
}

svec2_t get_uniform_vec2(uint32_t idx)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

    return(gl.uniform[idx].vec2);
}

svec3_t get_uniform_vec3(uint32_t idx)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

    return(gl.uniform[idx].vec3);
}

svec4_t get_uniform_vec4(uint32_t idx)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(idx < MAX_UNIFORMS);

    return(gl.uniform[idx].vec4);
}

void set_framebuffer(framebuffer_t *framebuffer)
{
    ASSERT(gl.valid); // call setup_gl() first!

    gl.framebuffer = framebuffer;
}

void set_fragment_shader(fragment_shader_t *shader)
{
    ASSERT(gl.valid); // call setup_gl() first!

    gl.fragment_shader = shader;
}

void clear_depth(framebuffer_t* framebuffer, float depth)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(framebuffer);

    uint32_t num = framebuffer->w*framebuffer->h;
    __stosd((unsigned long *)framebuffer->depth, *(unsigned long *)&depth, num);
}

void clear_color(framebuffer_t *framebuffer, vec4_t color)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(framebuffer);

    uint32_t color32 = pack_color(color);
    uint32_t num = framebuffer->w*framebuffer->h;
    __stosd((unsigned long *)framebuffer->color, *(unsigned long *)&color32, num);
}

static void apply_address_mode(texture_t *texture, int32_t *t, int32_t *s)
{
    switch(texture->address_mode)
    {
        case ADDRESS_MODE_REPEAT:
        {
            *t = *t % texture->w;
            *s = *s % texture->h;
            *t = (*t < 0) ? ((texture->w - 1) + *t) : *t;
            *s = (*s < 0) ? ((texture->h - 1) + *s) : *s;
            break;
        }
        case ADDRESS_MODE_CLAMP:
        {
            *t = CLAMP(*t, 0, (texture->w - 1));
            *s = CLAMP(*s, 0, (texture->h - 1));
            break;
        }
        default: ASSERT(!"Not implemented!"); break;
    }
}

svec4_t sample_texture(texture_t *texture, svec2_t uv)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(texture);

#if defined(FAST_MODE)
    int32_t ALIGNED(16) t[4];
    int32_t ALIGNED(16) s[4];
    int32_t ALIGNED(16) texel[4];

    _mm_store_si128((__m128i *)t, _mm_cvttps_epi32(_mm_mul_ps(uv.x, _mm_set1_ps((float)texture->w))));
    _mm_store_si128((__m128i *)s, _mm_cvttps_epi32(_mm_mul_ps(uv.y, _mm_set1_ps((float)texture->h))));

    apply_address_mode(texture, &t[0], &s[0]);
    apply_address_mode(texture, &t[1], &s[1]);
    apply_address_mode(texture, &t[2], &s[2]);
    apply_address_mode(texture, &t[3], &s[3]);

    texel[0] = s[0]*texture->w + t[0];
    texel[1] = s[1]*texture->w + t[1];
    texel[2] = s[2]*texture->w + t[2];
    texel[3] = s[3]*texture->w + t[3];

    uint32_t ALIGNED(16) color32[4];
    color32[0] = texture->pixels[texel[0]];
    color32[1] = texture->pixels[texel[1]];
    color32[2] = texture->pixels[texel[2]];
    color32[3] = texture->pixels[texel[3]];

    xvec4_t color = xunpack_color(_mm_load_si128((__m128i *)color32));
#else
    int32_t t = (int32_t)(uv.x*(float)texture->w);
    int32_t s = (int32_t)(uv.y*(float)texture->h);
    apply_address_mode(texture, &t, &s);
    int32_t texel = s*texture->w + t;
    vec4_t color = unpack_color(texture->pixels[texel]);
#endif

    return(color);
}

void draw(vertex_t *vertices, uint32_t *indices, uint32_t num_indices)
{
    ASSERT(gl.valid); // call setup_gl() first!
    ASSERT(gl.framebuffer); // call set_framebuffer() with a valid framebuffer first!
    ASSERT(gl.fragment_shader); // call set_fragment_shader() with a valid fragment shader first!

    gl.num_triangles = 0;

    vec2_t screen_size = { (float)gl.framebuffer->w, (float)gl.framebuffer->h };

    vec2_t screen_min = screen_size;
    vec2_t screen_max = { 0, 0 };

    m4x4_t mvp = m4x4_mul(m4x4_mul(gl.model, gl.view), gl.proj);

    for(uint32_t i = 0; i < num_indices; i += 3)
    {
        triangle_t triangles[2];
        uint32_t num_triangles = 0;

        {
            homogeneous_vertex_t v[3];
            convert_vertex(&vertices[indices[i]], &v[0]);
            convert_vertex(&vertices[indices[i + 1]], &v[1]);
            convert_vertex(&vertices[indices[i + 2]], &v[2]);

            v[0].hpos = m4x4_mulv(mvp, (vec4_t){ v[0].pos.x, v[0].pos.y, v[0].pos.z, 1.0f });
            v[1].hpos = m4x4_mulv(mvp, (vec4_t){ v[1].pos.x, v[1].pos.y, v[1].pos.z, 1.0f });
            v[2].hpos = m4x4_mulv(mvp, (vec4_t){ v[2].pos.x, v[2].pos.y, v[2].pos.z, 1.0f });

            v[0].normal = m4x4_mulv(gl.model, (vec4_t){ v[0].normal.x, v[0].normal.y, v[0].normal.z, 0.0f }).xyz;
            v[1].normal = m4x4_mulv(gl.model, (vec4_t){ v[1].normal.x, v[1].normal.y, v[1].normal.z, 0.0f }).xyz;
            v[2].normal = m4x4_mulv(gl.model, (vec4_t){ v[2].normal.x, v[2].normal.y, v[2].normal.z, 0.0f }).xyz;

            v[0].pos = m4x4_mulv(gl.model, (vec4_t){ v[0].pos.x, v[0].pos.y, v[0].pos.z, 1.0f }).xyz;
            v[1].pos = m4x4_mulv(gl.model, (vec4_t){ v[1].pos.x, v[1].pos.y, v[1].pos.z, 1.0f }).xyz;
            v[2].pos = m4x4_mulv(gl.model, (vec4_t){ v[2].pos.x, v[2].pos.y, v[2].pos.z, 1.0f }).xyz;

            // TODO: We currently only clip against the near plane. This might not be robust enough in all cases.
            uint8_t clip = 0;
            clip |= (uint8_t)(v[0].hpos.w <= 0.0f || (v[0].hpos.z < -v[0].hpos.w)) << 0;
            clip |= (uint8_t)(v[1].hpos.w <= 0.0f || (v[1].hpos.z < -v[1].hpos.w)) << 1;
            clip |= (uint8_t)(v[2].hpos.w <= 0.0f || (v[2].hpos.z < -v[2].hpos.w)) << 2;

            switch(clip)
            {
                case 0b000:
                {
                    memcpy(&triangles[num_triangles++], v, sizeof(triangle_t));
                    break;
                }
                case 0b001:
                {
                    homogeneous_vertex_t v01 = clip_edge(&v[1], &v[0]);
                    homogeneous_vertex_t v02 = clip_edge(&v[2], &v[0]);
                    triangles[num_triangles++] = (triangle_t){ v01, v[1], v[2] };
                    triangles[num_triangles++] = (triangle_t){ v01, v[2], v02 };
                    break;
                }
                case 0b010:
                {
                    homogeneous_vertex_t v12 = clip_edge(&v[2], &v[1]);
                    homogeneous_vertex_t v10 = clip_edge(&v[0], &v[1]);
                    triangles[num_triangles++] = (triangle_t){ v[0], v10, v[2] };
                    triangles[num_triangles++] = (triangle_t){ v10, v12, v[2] };
                    break;
                }
                case 0b100:
                {
                    homogeneous_vertex_t v20 = clip_edge(&v[0], &v[2]);
                    homogeneous_vertex_t v21 = clip_edge(&v[1], &v[2]);
                    triangles[num_triangles++] = (triangle_t){ v[0], v[1], v20 };
                    triangles[num_triangles++] = (triangle_t){ v20, v[1], v21 };
                    break;
                }
                case 0b011:
                {
                    v[0] = clip_edge(&v[2], &v[0]);
                    v[1] = clip_edge(&v[2], &v[1]);
                    memcpy(&triangles[num_triangles++], v, sizeof(triangle_t));
                    break;
                }
                case 0b101:
                {
                    v[0] = clip_edge(&v[1], &v[0]);
                    v[2] = clip_edge(&v[1], &v[2]);
                    memcpy(&triangles[num_triangles++], v, sizeof(triangle_t));
                    break;
                }
                case 0b110:
                {
                    v[1] = clip_edge(&v[0], &v[1]);
                    v[2] = clip_edge(&v[0], &v[2]);
                    memcpy(&triangles[num_triangles++], v, sizeof(triangle_t));
                    break;
                }
            }
        }

        for(uint32_t j = 0; j < num_triangles; j++)
        {
            triangle_t *triangle = &triangles[j];

            homogeneous_vertex_t *v0 = &triangle->v0;
            homogeneous_vertex_t *v1 = &triangle->v1;
            homogeneous_vertex_t *v2 = &triangle->v2;

            v0->hpos.w = 1.0f/v0->hpos.w;
            v1->hpos.w = 1.0f/v1->hpos.w;
            v2->hpos.w = 1.0f/v2->hpos.w;

            v0->hpos.xyz = vec3_muls(v0->hpos.xyz, v0->hpos.w);
            v1->hpos.xyz = vec3_muls(v1->hpos.xyz, v1->hpos.w);
            v2->hpos.xyz = vec3_muls(v2->hpos.xyz, v2->hpos.w);

            v0->hpos.x = ((v0->hpos.x*0.5f) + 0.5f)*screen_size.x;
            v1->hpos.x = ((v1->hpos.x*0.5f) + 0.5f)*screen_size.x;
            v2->hpos.x = ((v2->hpos.x*0.5f) + 0.5f)*screen_size.x;

            v0->hpos.y = ((v0->hpos.y*-0.5f) + 0.5f)*screen_size.y;
            v1->hpos.y = ((v1->hpos.y*-0.5f) + 0.5f)*screen_size.y;
            v2->hpos.y = ((v2->hpos.y*-0.5f) + 0.5f)*screen_size.y;

            float area = edge_func(v0->hpos.xy, v1->hpos.xy, v2->hpos.xy);

            if(area > 0)
                continue;

            vec2_t bounds_min = vec2_max(vec2_min(vec2_min(v0->hpos.xy, v1->hpos.xy), v2->hpos.xy), (vec2_t){ 0, 0 });
            vec2_t bounds_max = vec2_min(vec2_max(vec2_max(v0->hpos.xy, v1->hpos.xy), v2->hpos.xy), screen_size);

            if(bounds_min.x > bounds_max.x || bounds_min.y > bounds_max.y)
                continue;

            screen_min = vec2_min(screen_min, bounds_min);
            screen_max = vec2_max(screen_max, bounds_max);

            ASSERT(gl.num_triangles < MAX_TRIANGLES);
            gl.triangles[gl.num_triangles++] = *triangle;
        }
    }

    if(gl.num_triangles > 0)
    {
        int32_t start_x = (int32_t)(floor(screen_min.x/16.0f)*16.0f);
        int32_t start_y = (int32_t)(floor(screen_min.y/16.0f)*16.0f);
        int32_t end_x = (int32_t)(ceilf(screen_max.x/16.0f)*16.0f);
        int32_t end_y = (int32_t)(ceilf(screen_max.y/16.0f)*16.0f);

        // TODO: handle screen sizes that are not a multiple of four.
        int32_t dx = (end_x - start_x)/4;
        int32_t dy = (end_y - start_y)/4;

        int32_t num = 0;
        rect_t regions[16];
        for(int32_t x = start_x; x < end_x; x += dx)
        {
            for(int32_t y = start_y; y < end_y; y += dy)
            {
                int32_t w = MAX(MIN((x + dx), gl.framebuffer->w) - x, 0);
                int32_t h = MAX(MIN((y + dy), gl.framebuffer->h) - y, 0);
                if(w > 0 && h > 0)
                {
                    int32_t i = num++;
                    regions[i] = (rect_t){ x, y, w, h };
                    add_job(&draw_region, &regions[i]);
                }
            }
        }

        wait_for_all_jobs();
    }
}

void setup_gl(void)
{
    if(!gl.valid)
    {
        gl.valid = true;

        set_model_matrix(m4x4_identity());
        set_view_matrix(m4x4_identity());
        set_projection_matrix(m4x4_identity());

        setup_job_queue();
    }
}
