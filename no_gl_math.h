#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#define HALF_PI 1.57079632679f
#define PI 3.14159265359f
#define TAU 6.28318530718f
#define EPS 1.0e-6f

typedef union vec2_t
{
    struct
    {
        float x, y;
    };
    struct
    {
        float r, g;
    };
    float e[2];
} vec2_t;

typedef union xvec2_t
{
    struct
    {
        __m128 x, y;
    };
    struct
    {
        __m128 r, g;
    };    
    __m128 e[2];
} xvec2_t;

typedef union vec3_t
{
    struct
    {
        float x, y, z;
    };
    struct
    {
        vec2_t xy;
        PAD(4);
    };
    struct
    {
        PAD(4);
        vec2_t yz;
    };
    struct
    {
        float r, g, b;
    };
    struct
    {
        vec2_t rg;
        PAD(4);
    };
    struct
    {
        PAD(4);
        vec2_t gb;
    };
    float e[3];
} vec3_t;

typedef union xvec3_t
{
    struct
    {
        __m128 x, y, z;
    };
    struct
    {
        xvec2_t xy;
        PAD(16);
    };
    struct
    {
        PAD(16);
        xvec2_t yz;
    };
    struct
    {
        __m128 r, g, b;
    };
    struct
    {
        xvec2_t rg;
        PAD(16);
    };
    struct
    {
        PAD(16);
        xvec2_t gb;
    };
    __m128 e[3];
} xvec3_t;

typedef union vec4_t
{
    struct
    {
        float x, y, z, w;
    };
    struct
    {
        vec2_t xy, zw;
    };
    struct
    {
        PAD(4);
        vec2_t yz;
        PAD(4);
    };
    struct
    {
        vec3_t xyz;
        PAD(4);
    };
    struct
    {
        PAD(4);
        vec3_t yzw;
    };
    struct
    {
        float r, g, b, a;
    };
    struct
    {
        vec2_t rg, ba;
    };
    struct
    {
        PAD(4);
        vec2_t gb;
        PAD(4);
    };
    struct
    {
        vec3_t rgb;
        PAD(4);
    };
    struct
    {
        PAD(4);
        vec3_t gba;
    };
    float e[4];
} vec4_t;

typedef union xvec4_t
{
    struct
    {
        __m128 x, y, z, w;
    };
    struct
    {
        xvec2_t xy, zw;
    };
    struct
    {
        PAD(16);
        xvec2_t yz;
        PAD(16);
    };
    struct
    {
        xvec3_t xyz;
        PAD(16);
    };
    struct
    {
        PAD(16);
        xvec3_t yzw;
    };
    struct
    {
        __m128 r, g, b, a;
    };
    struct
    {
        xvec2_t rg, ba;
    };
    struct
    {
        PAD(16);
        xvec2_t gb;
        PAD(16);
    };
    struct
    {
        xvec3_t rgb;
        PAD(16);
    };
    struct
    {
        PAD(16);
        xvec3_t gba;
    };
    __m128 e[4];
} xvec4_t;

typedef struct m3x3_t
{
    float e[3][3];
} m3x3_t;

typedef union xm3x3_t
{
    __m128 e[3][3];
} xm3x3_t;

typedef struct m4x4_t
{
    float e[4][4];
} m4x4_t;

typedef union xm4x4_t
{
    __m128 e[4][4];
} xm4x4_t;

#if defined(FAST_MODE)

typedef __m128 scalar_t;
typedef xvec2_t svec2_t;
typedef xvec3_t svec3_t;
typedef xvec4_t svec4_t;
typedef xm3x3_t sm3x3_t;
typedef xm4x4_t sm4x4_t;

#define scalar(x) _mm_set_ps1(x)

#define smin _mm_min_ps
#define smax _mm_max_ps
#define spow xfast_pow

#define svec2_set1 xvec2_set1
#define svec2_set xvec2_set
#define svec2_normalize xvec2_normalize
#define svec2_sub xvec2_sub
#define svec2_add xvec2_add
#define svec2_mul xvec2_mul
#define svec2_muls xvec2_muls
#define svec2_divs xvec2_divs
#define svec2_length_sqr xvec2_length_sqr
#define svec2_length xvec2_length
#define svec2_normalize xvec2_normalize
#define svec2_lerp xvec2_lerp
#define svec2_pow xvec2_pow
#define svec2_dot xvec2_dot
#define svec2_cross xvec2_cross
#define svec2_min xvec2_min
#define svec2_max xvec2_max

#define svec3_set1 xvec3_set1
#define svec3_set xvec3_set
#define svec3_normalize xvec3_normalize
#define svec3_sub xvec3_sub
#define svec3_add xvec3_add
#define svec3_mul xvec3_mul
#define svec3_muls xvec3_muls
#define svec3_divs xvec3_divs
#define svec3_length_sqr xvec3_length_sqr
#define svec3_length xvec3_length
#define svec3_normalize xvec3_normalize
#define svec3_lerp xvec3_lerp
#define svec3_pow xvec3_pow
#define svec3_dot xvec3_dot
#define svec3_cross xvec3_cross
#define svec3_min xvec3_min
#define svec3_max xvec3_max

#define svec4_set1 xvec4_set1
#define svec4_set xvec4_set
#define svec4_normalize xvec4_normalize
#define svec4_sub xvec4_sub
#define svec4_add xvec4_add
#define svec4_mul xvec4_mul
#define svec4_muls xvec4_muls
#define svec4_divs xvec4_divs
#define svec4_length_sqr xvec4_length_sqr
#define svec4_length xvec4_length
#define svec4_normalize xvec4_normalize
#define svec4_lerp xvec4_lerp
#define svec4_pow xvec4_pow
#define svec4_dot xvec4_dot
#define svec4_cross xvec4_cross
#define svec4_min xvec4_min
#define svec4_max xvec4_max
    
#define sm3x3_mulv xm3x3_mulv
    
#define sm4x4_mulv xm4x4_mulv

#else

typedef float scalar_t;
typedef vec2_t svec2_t;
typedef vec3_t svec3_t;
typedef vec4_t svec4_t;
typedef m3x3_t sm3x3_t;
typedef m4x4_t sm4x4_t;

#define scalar(x) x

#define smin MIN
#define smax MAX
#define spow fast_pow

#define svec2_set1 vec2_set1
#define svec2_set vec2_set
#define svec2_normalize vec2_normalize
#define svec2_sub vec2_sub
#define svec2_add vec2_add
#define svec2_mul vec2_mul
#define svec2_muls vec2_muls
#define svec2_divs vec2_divs
#define svec2_length_sqr vec2_length_sqr
#define svec2_length vec2_length
#define svec2_normalize vec2_normalize
#define svec2_lerp vec2_lerp
#define svec2_pow vec2_pow
#define svec2_dot vec2_dot
#define svec2_cross vec2_cross
#define svec2_min vec2_min
#define svec2_max vec2_max

#define svec3_set1 vec3_set1
#define svec3_set vec3_set
#define svec3_normalize vec3_normalize
#define svec3_sub vec3_sub
#define svec3_add vec3_add
#define svec3_mul vec3_mul
#define svec3_muls vec3_muls
#define svec3_divs vec3_divs
#define svec3_length_sqr vec3_length_sqr
#define svec3_length vec3_length
#define svec3_normalize vec3_normalize
#define svec3_lerp vec3_lerp
#define svec3_pow vec3_pow
#define svec3_dot vec3_dot
#define svec3_cross vec3_cross
#define svec3_min vec3_min
#define svec3_max vec3_max

#define svec4_set1 vec4_set1
#define svec4_set vec4_set
#define svec4_normalize vec4_normalize
#define svec4_sub vec4_sub
#define svec4_add vec4_add
#define svec4_mul vec4_mul
#define svec4_muls vec4_muls
#define svec4_divs vec4_divs
#define svec4_length_sqr vec4_length_sqr
#define svec4_length vec4_length
#define svec4_normalize vec4_normalize
#define svec4_lerp vec4_lerp
#define svec4_pow vec4_pow
#define svec4_dot vec4_dot
#define svec4_cross vec4_cross
#define svec4_min vec4_min
#define svec4_max vec4_max

#define sm3x3_mulv m3x3_mulv

#define sm4x4_mulv m4x4_mulv

#endif

static float clamp(float v, float a, float b)
{
    return(CLAMP(v, a, b));
}

static float clamp01(float v)
{
    return(clamp(v, 0, 1));
}

static float lerp(float a, float b, float t)
{
    float result = a*(1-t)+b*t;
    return(result);
}

static float deg_to_rad(float angle)
{
    float result = (angle/180.0f)*PI;
    return(result);
}

static float rad_to_deg(float angle)
{
    float result = (angle/PI)*180.0f;
    return(result);
}

// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
static float fast_pow(float value, float power)
{
    double a = (double)value;
    double b = (double)power;
    // calculate approximation with fraction of the exponent
    int e = (int) b;
    union
    {
        double d;
        int x[2];
    } u = { a };
    u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;

    // exponentiation by squaring with the exponent's integer part
    // double r = u.d makes everything much slower, not sure why
    double r = 1.0;
    while (e)
    {
        if (e & 1)
        {
            r *= a;
        }
        a *= a;
        e >>= 1;
    }
    return((float)(r * u.d));
}

static uint64_t next_pow2(uint64_t value)
{
    uint64_t pow2 = value;
    pow2--;
    pow2 |= pow2 >> 1;
    pow2 |= pow2 >> 2;
    pow2 |= pow2 >> 4;
    pow2 |= pow2 >> 8;
    pow2 |= pow2 >> 16;
    pow2++;
    return(pow2);
}

static vec2_t vec2_set1(float v)
{
    vec2_t result;
    result.x = v;
    result.y = v;
    return(result);
}

static vec2_t vec2_set(float x, float y)
{
    vec2_t result;
    result.x = x;
    result.y = y;
    return(result);
}

static vec2_t vec2_add(vec2_t a, vec2_t b)
{
    vec2_t result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    return(result);
}

static vec2_t vec2_sub(vec2_t a, vec2_t b)
{
    vec2_t result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    return(result);
}

static vec2_t vec2_mul(vec2_t a, vec2_t b)
{
    vec2_t result;
    result.x = a.x*b.x;
    result.y = a.y*b.y;
    return(result);
}

static vec2_t vec2_muls(vec2_t a, float b)
{
    vec2_t result;
    result.x = a.x*b;
    result.y = a.y*b;
    return(result);
}

static vec2_t vec2_divs(vec2_t a, float b)
{
    vec2_t result;
    result.x = a.x/b;
    result.y = a.y/b;
    return(result);
}

static vec2_t vec2_abs(vec2_t v)
{
    vec2_t result;
    result.x = fabsf(v.x);
    result.y = fabsf(v.y);
    return(result);
}

static bool vec2_equal(vec2_t a, vec2_t b)
{
    bool result = ((a.x == b.x) && (a.y == b.y));
    return(result);
}

static float vec2_length_sqr(vec2_t v)
{
    float result = (v.x * v.x) + (v.y * v.y);
    return(result);
}

static float vec2_length(vec2_t v)
{
    float result = sqrtf(vec2_length_sqr(v));
    return(result);
}

static vec2_t vec2_normalize(vec2_t v)
{
    vec2_t result = vec2_divs(v, vec2_length(v));
    return(result);
}

static vec2_t vec2_lerp(vec2_t a, vec2_t b, float t)
{
    vec2_t result = vec2_add(vec2_muls(a, 1.0f - t), vec2_muls(b, t));
    return(result);
}

static vec2_t vec2_pow(vec2_t a, float b)
{
    vec2_t result;
    result.x = fast_pow(a.x, b);
    result.y = fast_pow(a.y, b);
    return(result);
}

static float vec2_dot(vec2_t a, vec2_t b)
{
    float result = (a.x * b.x) + (a.y * b.y);
    return(result);
}

static float vec2_cross(vec2_t a, vec2_t b)
{
    float result = (a.x * b.y) - (b.x * a.y);
    return(result);
}

static vec2_t vec2_min(vec2_t a, vec2_t b)
{
    vec2_t result;
    result.x = MIN(a.x, b.x);
    result.y = MIN(a.y, b.y);
    return(result);
}

static vec2_t vec2_max(vec2_t a, vec2_t b)
{
    vec2_t result;
    result.x = MAX(a.x, b.x);
    result.y = MAX(a.y, b.y);
    return(result);
}

static vec3_t vec3_set1(float v)
{
    vec3_t result;
    result.x = v;
    result.y = v;
    result.z = v;
    return(result);
}

static vec3_t vec3_set(float x, float y, float z)
{
    vec3_t result;
    result.x = x;
    result.y = y;
    result.z = z;
    return(result);
}

static vec3_t vec3_add(vec3_t a, vec3_t b)
{
    vec3_t result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return(result);
}

static vec3_t vec3_sub(vec3_t a, vec3_t b)
{
    vec3_t result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return(result);
}

static vec3_t vec3_mul(vec3_t a, vec3_t b)
{
    vec3_t result;
    result.x = a.x * b.x;
    result.y = a.y * b.y;
    result.z = a.z * b.z;
    return(result);
}

static vec3_t vec3_muls(vec3_t a, float b)
{
    vec3_t result;
    result.x = a.x * b;
    result.y = a.y * b;
    result.z = a.z * b;
    return(result);
}

static vec3_t vec3_divs(vec3_t a, float b)
{
    vec3_t result;
    result.x = a.x / b;
    result.y = a.y / b;
    result.z = a.z / b;
    return(result);
}

static vec3_t vec3_abs(vec3_t v)
{
    vec3_t result;
    result.x = fabsf(v.x);
    result.y = fabsf(v.y);
    result.z = fabsf(v.z);
    return(result);
}

static bool vec3_equal(vec3_t a, vec3_t b)
{
    bool result = ((a.x == b.x) && (a.y == b.y) && (a.z == b.z));
    return(result);
}

static float vec3_length_sqr(vec3_t v)
{
    float result = (v.x * v.x) + (v.y * v.y) + (v.z * v.z);
    return(result);
}

static float vec3_length(vec3_t v)
{
    float result = sqrtf(vec3_length_sqr(v));
    return(result);
}

static vec3_t vec3_normalize(vec3_t v)
{
    vec3_t result = vec3_divs(v, vec3_length(v));
    return(result);
}

static vec3_t vec3_lerp(vec3_t a, vec3_t b, float t)
{
    vec3_t result = vec3_add(vec3_muls(a, 1.0f - t), vec3_muls(b, t));
    return(result);
}

static vec3_t vec3_pow(vec3_t a, float b)
{
    vec3_t result;
    result.x = fast_pow(a.x, b);
    result.y = fast_pow(a.y, b);
    result.z = fast_pow(a.z, b);
    return(result);
}

static float vec3_dot(vec3_t a, vec3_t b)
{
    float result = (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    return(result);
}

static vec3_t vec3_cross(vec3_t a, vec3_t b)
{
    vec3_t result;
    result.x = (a.y * b.z) - (b.y * a.z);
    result.y = (a.z * b.x) - (b.z * a.x);
    result.z = (a.x * b.y) - (b.x * a.y);
    return(result);
}

static vec3_t vec3_min(vec3_t a, vec3_t b)
{
    vec3_t result;
    result.x = MIN(a.x, b.x);
    result.y = MIN(a.y, b.y);
    result.z = MIN(a.z, b.z);
    return(result);
}

static vec3_t vec3_max(vec3_t a, vec3_t b)
{
    vec3_t result;
    result.x = MAX(a.x, b.x);
    result.y = MAX(a.y, b.y);
    result.z = MAX(a.z, b.z);
    return(result);
}

static vec4_t vec4_set1(float v)
{
    vec4_t result;
    result.x = v;
    result.y = v;
    result.z = v;
    result.w = v;
    return(result);
}

static vec4_t vec4_set(float x, float y, float z, float w)
{
    vec4_t result;
    result.x = x;
    result.y = y;
    result.z = z;
    result.w = w;
    return(result);
}

static vec4_t vec4_add(vec4_t a, vec4_t b)
{
    vec4_t result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    result.w = a.w + b.w;
    return(result);
}

static vec4_t vec4_sub(vec4_t a, vec4_t b)
{
    vec4_t result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    result.w = a.w - b.w;
    return(result);
}

static vec4_t vec4_mul(vec4_t a, vec4_t b)
{
    vec4_t result;
    result.x = a.x * b.x;
    result.y = a.y * b.y;
    result.z = a.z * b.z;
    result.w = a.w * b.w;
    return(result);
}

static vec4_t vec4_muls(vec4_t a, float b)
{
    vec4_t result;
    result.x = a.x * b;
    result.y = a.y * b;
    result.z = a.z * b;
    result.w = a.w * b;
    return(result);
}

static vec4_t vec4_divs(vec4_t a, float b)
{
    vec4_t result;
    result.x = a.x / b;
    result.y = a.y / b;
    result.z = a.z / b;
    result.w = a.w / b;
    return(result);
}

static vec4_t vec4_abs(vec4_t v)
{
    vec4_t result;
    result.x = fabsf(v.x);
    result.y = fabsf(v.y);
    result.z = fabsf(v.z);
    result.w = fabsf(v.w);
    return(result);
}

static bool vec4_equal(vec4_t a, vec4_t b)
{
    bool result = ((a.x == b.x) && (a.y == b.y) && (a.z == b.z) && (a.w == b.w));
    return(result);
}

static float vec4_dot(vec4_t a, vec4_t b)
{
    float result = (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
    return(result);
}

static float vec4_length_sqr(vec4_t a)
{
    float result = vec4_dot(a, a);
    return(result);
}

static float vec4_length(vec4_t a)
{
    float result = sqrtf(vec4_dot(a, a));
    return(result);
}

static vec4_t vec4_normalize(vec4_t a)
{
    vec4_t result = vec4_divs(a, vec4_length(a));
    return(result);
}

static vec4_t vec4_lerp(vec4_t a, vec4_t b, float t)
{
    vec4_t result = vec4_add(vec4_muls(a, 1.0f - t), vec4_muls(b, t));
    return(result);
}

static vec4_t vec4_pow(vec4_t a, float b)
{
    vec4_t result;
    result.x = fast_pow(a.x, b);
    result.y = fast_pow(a.y, b);
    result.z = fast_pow(a.z, b);
    result.w = fast_pow(a.w, b);
    return(result);
}

static vec4_t vec4_cross(vec4_t a, vec4_t b)
{
    vec4_t result;
    result.x = (a.y * b.z) - (b.y * a.z);
    result.y = (a.z * b.x) - (b.z * a.x);
    result.z = (a.x * b.y) - (b.x * a.y);
    result.w = 0;
    return(result);
}

static vec4_t vec4_min(vec4_t a, vec4_t b)
{
    vec4_t result;
    result.x = MIN(a.x, b.x);
    result.y = MIN(a.y, b.y);
    result.z = MIN(a.z, b.z);
    result.w = MIN(a.w, b.w);
    return(result);
}

static vec4_t vec4_max(vec4_t a, vec4_t b)
{
    vec4_t result;
    result.x = MAX(a.x, b.x);
    result.y = MAX(a.y, b.y);
    result.z = MAX(a.z, b.z);
    result.w = MAX(a.w, b.w);
    return(result);
}

static vec4_t vec4_reflect(vec4_t n, vec4_t i)
{
    vec4_t result = vec4_sub(i, vec4_muls(n, 2.0f*vec4_dot(n, i)));
    return(result);
}

static m3x3_t m3x3_identity(void)
{
    m3x3_t result =
        {
            {{1, 0, 0},
             {0, 1, 0},
             {0, 0, 1}}
        };
    return(result);
}

static m3x3_t m3x3_transpose(m3x3_t m)
{
    m3x3_t result;
    for(uint32_t i = 0; i < 3; i++)
        for(uint32_t j = 0; j < 3; j++)
            result.e[i][j] = m.e[j][i];

    return(result);
}

static m3x3_t m3x3_mul(m3x3_t a, m3x3_t b)
{
    m3x3_t result;

    result.e[0][0] = a.e[0][0]*b.e[0][0] + a.e[0][1]*b.e[1][0] + a.e[0][2]*b.e[2][0];
    result.e[0][1] = a.e[0][0]*b.e[0][1] + a.e[0][1]*b.e[1][1] + a.e[0][2]*b.e[2][1];
    result.e[0][2] = a.e[0][0]*b.e[0][2] + a.e[0][1]*b.e[1][2] + a.e[0][2]*b.e[2][2];

    result.e[1][0] = a.e[1][0]*b.e[0][0] + a.e[1][1]*b.e[1][0] + a.e[1][2]*b.e[2][0];
    result.e[1][1] = a.e[1][0]*b.e[0][1] + a.e[1][1]*b.e[1][1] + a.e[1][2]*b.e[2][1];
    result.e[1][2] = a.e[1][0]*b.e[0][2] + a.e[1][1]*b.e[1][2] + a.e[1][2]*b.e[2][2];

    result.e[2][0] = a.e[2][0]*b.e[0][0] + a.e[2][1]*b.e[1][0] + a.e[2][2]*b.e[2][0];
    result.e[2][1] = a.e[2][0]*b.e[0][1] + a.e[2][1]*b.e[1][1] + a.e[2][2]*b.e[2][1];
    result.e[2][2] = a.e[2][0]*b.e[0][2] + a.e[2][1]*b.e[1][2] + a.e[2][2]*b.e[2][2];

    return(result);
}

static vec3_t m3x3_mulv(m3x3_t a, vec3_t b)
{
    vec3_t result;
    result.x = a.e[0][0]*b.x + a.e[1][0]*b.y + a.e[2][0]*b.z;
    result.y = a.e[0][1]*b.x + a.e[1][1]*b.y + a.e[2][1]*b.z;
    result.z = a.e[0][2]*b.x + a.e[1][2]*b.y + a.e[2][2]*b.z;
    return(result);
}

static m3x3_t m3x3_inverse(m3x3_t m)
{
    m3x3_t result = m;

    float e00 = (m.e[1][1]*m.e[2][2] - m.e[1][2]*m.e[2][1]);
    float e01 = -(m.e[1][0]*m.e[2][2] - m.e[1][2]*m.e[2][0]);
    float e02 = (m.e[1][0]*m.e[2][1] - m.e[1][1]*m.e[2][0]);
    float det = m.e[0][0]*e00 + m.e[0][1]*e01 + m.e[0][2]*e02;

    if(det != 0)
    {
        det = 1/det;

        result.e[0][0] = det*e00;
        result.e[1][0] = det*e01;
        result.e[2][0] = det*e02;

        result.e[0][1] = -det*(m.e[0][1]*m.e[2][2] - m.e[0][2]*m.e[2][1]);
        result.e[1][1] = det*(m.e[0][0]*m.e[2][2] - m.e[0][2]*m.e[2][0]);
        result.e[2][1] = -det*(m.e[0][0]*m.e[2][1] - m.e[0][1]*m.e[2][0]);

        result.e[0][2] = det*(m.e[0][1]*m.e[1][2] - m.e[0][2]*m.e[1][1]);
        result.e[1][2] = -det*(m.e[0][0]*m.e[1][2] - m.e[0][2]*m.e[1][0]);
        result.e[2][2] = det*(m.e[0][0]*m.e[1][1] - m.e[0][1]*m.e[1][0]);
    }

    return(result);
}

static m4x4_t m4x4_identity(void)
{
    m4x4_t result =
        {
            {{1, 0, 0, 0},
             {0, 1, 0, 0},
             {0, 0, 1, 0},
             {0, 0, 0, 1}}
        };
    return(result);
}

static m4x4_t m4x4_translate(vec3_t pos)
{
    m4x4_t result =
        {
            {{1, 0, 0, pos.x},
             {0, 1, 0, pos.y},
             {0, 0, 1, pos.z},
             {0, 0, 0, 1}}
        };
    return(result);
}

static m4x4_t m4x4_rotate(vec3_t rot)
{
    float cx = cosf(rot.x);
    float sx = sinf(rot.x);
    float cy = cosf(rot.y);
    float sy = sinf(rot.y);
    float cz = cosf(rot.z);
    float sz = sinf(rot.z);
    float sxsy = (sx * sy);
    float cxcy = (cx * sy);

    m4x4_t result;
    result.e[0][0] = (cy * cz);
    result.e[0][1] = (cy * sz);
    result.e[0][2] = (-sy);
    result.e[0][3] = 0;
    result.e[1][0] = ((sxsy * cz) - (cx * sz));
    result.e[1][1] = ((sxsy * sz) + (cx * cz));
    result.e[1][2] = (sx * cy);
    result.e[1][3] = 0;
    result.e[2][0] = ((cxcy * cz) + (sx * sz));
    result.e[2][1] = ((cxcy * sz) - (sx * cz));
    result.e[2][2] = (cx * cy);
    result.e[2][3] = 0;
    result.e[3][0] = 0;
    result.e[3][1] = 0;
    result.e[3][2] = 0;
    result.e[3][3] = 1.0f;

    return(result);
}

static m4x4_t m4x4_scale(vec3_t scale)
{
    m4x4_t mat =
        {
            {{scale.x, 0, 0, 0},
             {0, scale.y, 0, 0},
             {0, 0, scale.z, 0},
             {0, 0, 0, 1}}
        };
    return(mat);
}

static m4x4_t m4x4_look_at(vec3_t eye, vec3_t at, vec3_t up)
{
    vec3_t z_axis = vec3_normalize(vec3_sub(at, eye));
    vec3_t x_axis = vec3_normalize(vec3_cross(up, z_axis));
    vec3_t y_axis = vec3_cross(z_axis, x_axis);

    m4x4_t result;
    result.e[0][0] = x_axis.x;
    result.e[0][1] = y_axis.x;
    result.e[0][2] = z_axis.x;
    result.e[0][3] = 0;
    result.e[1][0] = x_axis.y;
    result.e[1][1] = y_axis.y;
    result.e[1][2] = z_axis.y;
    result.e[1][3] = 0;
    result.e[2][0] = x_axis.z;
    result.e[2][1] = y_axis.z;
    result.e[2][2] = z_axis.z;
    result.e[2][3] = 0;
    result.e[3][0] = -vec3_dot(x_axis, eye);
    result.e[3][1] = -vec3_dot(y_axis, eye);
    result.e[3][2] = -vec3_dot(z_axis, eye);
    result.e[3][3] = 1;

    return result;
}

static m4x4_t m4x4_perspective(float fov_y, float aspect, float nz, float fz)
{
    float h = (1.0f/tanf(fov_y * 0.5f));
    float w = (h/aspect);

    m4x4_t result;
    result.e[0][0] = w;
    result.e[0][1] = 0;
    result.e[0][2] = 0;
    result.e[0][3] = 0;
    result.e[1][0] = 0;
    result.e[1][1] = h;
    result.e[1][2] = 0;
    result.e[1][3] = 0;
    result.e[2][0] = 0;
    result.e[2][1] = 0;
    result.e[2][2] = (-(nz+fz)/(nz-fz));
    result.e[2][3] = 1.0f;
    result.e[3][0] = 0;
    result.e[3][1] = 0;
    result.e[3][2] = ((2.0f*fz*nz)/(nz-fz));
    result.e[3][3] = 0;

    return(result);
}

static m4x4_t m4x4_orthographic(float w, float h, float nz, float fz)
{
    m4x4_t result;
    result.e[0][0] = 2.0f/w;
    result.e[0][1] = 0;
    result.e[0][2] = 0;
    result.e[0][3] = 0;
    result.e[1][0] = 0;
    result.e[1][1] = -2.0f/h;
    result.e[1][2] = 0;
    result.e[1][3] = 0;
    result.e[2][0] = 0;
    result.e[2][1] = 0;
    result.e[2][2] = 1.0f/(fz-nz);
    result.e[2][3] = 0;
    result.e[3][0] = 0;
    result.e[3][1] = 0;
    result.e[3][2] = nz/(nz-fz);
    result.e[3][3] = 1;

    return result;
}

static m4x4_t m4x4_inverse(m4x4_t m)
{
    m4x4_t result = m;
    float det = (m.e[0][0] * m.e[1][1] - m.e[0][1] * m.e[1][0]) * (m.e[2][2] * m.e[3][3] - m.e[2][3] * m.e[3][2]) -
        (m.e[0][0] * m.e[1][2] - m.e[0][2] * m.e[1][0]) * (m.e[2][1] * m.e[3][3] - m.e[2][3] * m.e[3][1]) +
        (m.e[0][0] * m.e[1][3] - m.e[0][3] * m.e[1][0]) * (m.e[2][1] * m.e[3][2] - m.e[2][2] * m.e[3][1]) +
        (m.e[0][1] * m.e[1][2] - m.e[0][2] * m.e[1][1]) * (m.e[2][0] * m.e[3][3] - m.e[2][3] * m.e[3][0]) -
        (m.e[0][1] * m.e[1][3] - m.e[0][3] * m.e[1][1]) * (m.e[2][0] * m.e[3][2] - m.e[2][2] * m.e[3][0]) +
        (m.e[0][2] * m.e[1][3] - m.e[0][3] * m.e[1][2]) * (m.e[2][0] * m.e[3][1] - m.e[2][1] * m.e[3][0]);

    if(det != 0)
    {
        det = 1/det;

        result.e[0][0] = det * (m.e[1][1] * (m.e[2][2] * m.e[3][3] - m.e[2][3] * m.e[3][2]) +
                                m.e[1][2] * (m.e[2][3] * m.e[3][1] - m.e[2][1] * m.e[3][3]) +
                                m.e[1][3] * (m.e[2][1] * m.e[3][2] - m.e[2][2] * m.e[3][1]));
        result.e[0][1] = det * (m.e[2][1] * (m.e[0][2] * m.e[3][3] - m.e[0][3] * m.e[3][2]) +
                                m.e[2][2] * (m.e[0][3] * m.e[3][1] - m.e[0][1] * m.e[3][3]) +
                                m.e[2][3] * (m.e[0][1] * m.e[3][2] - m.e[0][2] * m.e[3][1]));
        result.e[0][2] = det * (m.e[3][1] * (m.e[0][2] * m.e[1][3] - m.e[0][3] * m.e[1][2]) +
                                m.e[3][2] * (m.e[0][3] * m.e[1][1] - m.e[0][1] * m.e[1][3]) +
                                m.e[3][3] * (m.e[0][1] * m.e[1][2] - m.e[0][2] * m.e[1][1]));
        result.e[0][3] = det * (m.e[0][1] * (m.e[1][3] * m.e[2][2] - m.e[1][2] * m.e[2][3]) +
                                m.e[0][2] * (m.e[1][1] * m.e[2][3] - m.e[1][3] * m.e[2][1]) +
                                m.e[0][3] * (m.e[1][2] * m.e[2][1] - m.e[1][1] * m.e[2][2]));
        result.e[1][0] = det * (m.e[1][2] * (m.e[2][0] * m.e[3][3] - m.e[2][3] * m.e[3][0]) +
                                m.e[1][3] * (m.e[2][2] * m.e[3][0] - m.e[2][0] * m.e[3][2]) +
                                m.e[1][0] * (m.e[2][3] * m.e[3][2] - m.e[2][2] * m.e[3][3]));
        result.e[1][1] = det * (m.e[2][2] * (m.e[0][0] * m.e[3][3] - m.e[0][3] * m.e[3][0]) +
                                m.e[2][3] * (m.e[0][2] * m.e[3][0] - m.e[0][0] * m.e[3][2]) +
                                m.e[2][0] * (m.e[0][3] * m.e[3][2] - m.e[0][2] * m.e[3][3]));
        result.e[1][2] = det * (m.e[3][2] * (m.e[0][0] * m.e[1][3] - m.e[0][3] * m.e[1][0]) +
                                m.e[3][3] * (m.e[0][2] * m.e[1][0] - m.e[0][0] * m.e[1][2]) +
                                m.e[3][0] * (m.e[0][3] * m.e[1][2] - m.e[0][2] * m.e[1][3]));
        result.e[1][3] = det * (m.e[0][2] * (m.e[1][3] * m.e[2][0] - m.e[1][0] * m.e[2][3]) +
                                m.e[0][3] * (m.e[1][0] * m.e[2][2] - m.e[1][2] * m.e[2][0]) +
                                m.e[0][0] * (m.e[1][2] * m.e[2][3] - m.e[1][3] * m.e[2][2]));
        result.e[2][0] = det * (m.e[1][3] * (m.e[2][0] * m.e[3][1] - m.e[2][1] * m.e[3][0]) +
                                m.e[1][0] * (m.e[2][1] * m.e[3][3] - m.e[2][3] * m.e[3][1]) +
                                m.e[1][1] * (m.e[2][3] * m.e[3][0] - m.e[2][0] * m.e[3][3]));
        result.e[2][1] = det * (m.e[2][3] * (m.e[0][0] * m.e[3][1] - m.e[0][1] * m.e[3][0]) +
                                m.e[2][0] * (m.e[0][1] * m.e[3][3] - m.e[0][3] * m.e[3][1]) +
                                m.e[2][1] * (m.e[0][3] * m.e[3][0] - m.e[0][0] * m.e[3][3]));
        result.e[2][2] = det * (m.e[3][3] * (m.e[0][0] * m.e[1][1] - m.e[0][1] * m.e[1][0]) +
                                m.e[3][0] * (m.e[0][1] * m.e[1][3] - m.e[0][3] * m.e[1][1]) +
                                m.e[3][1] * (m.e[0][3] * m.e[1][0] - m.e[0][0] * m.e[1][3]));
        result.e[2][3] = det * (m.e[0][3] * (m.e[1][1] * m.e[2][0] - m.e[1][0] * m.e[2][1]) +
                                m.e[0][0] * (m.e[1][3] * m.e[2][1] - m.e[1][1] * m.e[2][3]) +
                                m.e[0][1] * (m.e[1][0] * m.e[2][3] - m.e[1][3] * m.e[2][0]));
        result.e[3][0] = det * (m.e[1][0] * (m.e[2][2] * m.e[3][1] - m.e[2][1] * m.e[3][2]) +
                                m.e[1][1] * (m.e[2][0] * m.e[3][2] - m.e[2][2] * m.e[3][0]) +
                                m.e[1][2] * (m.e[2][1] * m.e[3][0] - m.e[2][0] * m.e[3][1]));
        result.e[3][1] = det * (m.e[2][0] * (m.e[0][2] * m.e[3][1] - m.e[0][1] * m.e[3][2]) +
                                m.e[2][1] * (m.e[0][0] * m.e[3][2] - m.e[0][2] * m.e[3][0]) +
                                m.e[2][2] * (m.e[0][1] * m.e[3][0] - m.e[0][0] * m.e[3][1]));
        result.e[3][2] = det * (m.e[3][0] * (m.e[0][2] * m.e[1][1] - m.e[0][1] * m.e[1][2]) +
                                m.e[3][1] * (m.e[0][0] * m.e[1][2] - m.e[0][2] * m.e[1][0]) +
                                m.e[3][2] * (m.e[0][1] * m.e[1][0] - m.e[0][0] * m.e[1][1]));
        result.e[3][3] = det * (m.e[0][0] * (m.e[1][1] * m.e[2][2] - m.e[1][2] * m.e[2][1]) +
                                m.e[0][1] * (m.e[1][2] * m.e[2][0] - m.e[1][0] * m.e[2][2]) +
                                m.e[0][2] * (m.e[1][0] * m.e[2][1] - m.e[1][1] * m.e[2][0]));
    }

    return(result);
}

static m4x4_t m4x4_mul(m4x4_t a, m4x4_t b)
{
    m4x4_t result;

    result.e[0][0] = (a.e[0][0]*b.e[0][0]) + (a.e[0][1]*b.e[1][0]) + (a.e[0][2]*b.e[2][0]) + (a.e[0][3]*b.e[3][0]);
    result.e[0][1] = (a.e[0][0]*b.e[0][1]) + (a.e[0][1]*b.e[1][1]) + (a.e[0][2]*b.e[2][1]) + (a.e[0][3]*b.e[3][1]);
    result.e[0][2] = (a.e[0][0]*b.e[0][2]) + (a.e[0][1]*b.e[1][2]) + (a.e[0][2]*b.e[2][2]) + (a.e[0][3]*b.e[3][2]);
    result.e[0][3] = (a.e[0][0]*b.e[0][3]) + (a.e[0][1]*b.e[1][3]) + (a.e[0][2]*b.e[2][3]) + (a.e[0][3]*b.e[3][3]);

    result.e[1][0] = (a.e[1][0]*b.e[0][0]) + (a.e[1][1]*b.e[1][0]) + (a.e[1][2]*b.e[2][0]) + (a.e[1][3]*b.e[3][0]);
    result.e[1][1] = (a.e[1][0]*b.e[0][1]) + (a.e[1][1]*b.e[1][1]) + (a.e[1][2]*b.e[2][1]) + (a.e[1][3]*b.e[3][1]);
    result.e[1][2] = (a.e[1][0]*b.e[0][2]) + (a.e[1][1]*b.e[1][2]) + (a.e[1][2]*b.e[2][2]) + (a.e[1][3]*b.e[3][2]);
    result.e[1][3] = (a.e[1][0]*b.e[0][3]) + (a.e[1][1]*b.e[1][3]) + (a.e[1][2]*b.e[2][3]) + (a.e[1][3]*b.e[3][3]);

    result.e[2][0] = (a.e[2][0]*b.e[0][0]) + (a.e[2][1]*b.e[1][0]) + (a.e[2][2]*b.e[2][0]) + (a.e[2][3]*b.e[3][0]);
    result.e[2][1] = (a.e[2][0]*b.e[0][1]) + (a.e[2][1]*b.e[1][1]) + (a.e[2][2]*b.e[2][1]) + (a.e[2][3]*b.e[3][1]);
    result.e[2][2] = (a.e[2][0]*b.e[0][2]) + (a.e[2][1]*b.e[1][2]) + (a.e[2][2]*b.e[2][2]) + (a.e[2][3]*b.e[3][2]);
    result.e[2][3] = (a.e[2][0]*b.e[0][3]) + (a.e[2][1]*b.e[1][3]) + (a.e[2][2]*b.e[2][3]) + (a.e[2][3]*b.e[3][3]);

    result.e[3][0] = (a.e[3][0]*b.e[0][0]) + (a.e[3][1]*b.e[1][0]) + (a.e[3][2]*b.e[2][0]) + (a.e[3][3]*b.e[3][0]);
    result.e[3][1] = (a.e[3][0]*b.e[0][1]) + (a.e[3][1]*b.e[1][1]) + (a.e[3][2]*b.e[2][1]) + (a.e[3][3]*b.e[3][1]);
    result.e[3][2] = (a.e[3][0]*b.e[0][2]) + (a.e[3][1]*b.e[1][2]) + (a.e[3][2]*b.e[2][2]) + (a.e[3][3]*b.e[3][2]);
    result.e[3][3] = (a.e[3][0]*b.e[0][3]) + (a.e[3][1]*b.e[1][3]) + (a.e[3][2]*b.e[2][3]) + (a.e[3][3]*b.e[3][3]);

    return(result);
}

static vec4_t m4x4_mulv(m4x4_t a, vec4_t b)
{
    vec4_t result;
    __m128 row0 = _mm_mul_ps(_mm_set_ps1(b.x), _mm_load_ps(&a.e[0][0]));
    __m128 row1 = _mm_mul_ps(_mm_set_ps1(b.y), _mm_load_ps(&a.e[1][0]));
    __m128 row2 = _mm_mul_ps(_mm_set_ps1(b.z), _mm_load_ps(&a.e[2][0]));
    __m128 row3 = _mm_mul_ps(_mm_set_ps1(b.w), _mm_load_ps(&a.e[3][0]));
    _mm_store_ps(&result.e[0], _mm_add_ps(_mm_add_ps(_mm_add_ps(row0, row1), row2), row3));
    return(result);
}

static m4x4_t m4x4_transpose(m4x4_t m)
{
    m4x4_t result;
    for(uint32_t i = 0; i < 4; i++)
        for(uint32_t j = 0; j < 4; j++)
            result.e[i][j] = m.e[j][i];

    return(result);
}

static __m128 xfast_pow(__m128 base, __m128 exponent)
{
    float ALIGNED(16) basev[4];
    float ALIGNED(16) exponentv[4];
    _mm_store_ps(basev, base);
    _mm_store_ps(exponentv, exponent);
    basev[0] = powf(basev[0], exponentv[0]);
    basev[1] = powf(basev[1], exponentv[1]);
    basev[2] = powf(basev[2], exponentv[2]);
    basev[3] = powf(basev[3], exponentv[3]);    
    return(_mm_load_ps(basev));
}

static xvec2_t xvec2_set1(__m128 v)
{
    xvec2_t result = { v, v };
    return(result);
}

static xvec2_t xvec2_set(__m128 x, __m128 y)
{
    xvec2_t result = { x, y };
    return(result);
}

static xvec2_t xvec2_setv(vec2_t v)
{
    xvec2_t result = { _mm_set_ps1(v.x), _mm_set_ps1(v.y) };
    return(result);
}

static xvec2_t xvec2_add(xvec2_t a, xvec2_t b)
{
    xvec2_t result;
    result.x = _mm_add_ps(a.x, b.x);
    result.y = _mm_add_ps(a.y, b.y);
    return(result);
}

static xvec2_t xvec2_sub(xvec2_t a, xvec2_t b)
{
    xvec2_t result;
    result.x = _mm_sub_ps(a.x, b.x);
    result.y = _mm_sub_ps(a.y, b.y);
    return(result);
}

static xvec2_t xvec2_mul(xvec2_t a, xvec2_t b)
{
    xvec2_t result;
    result.x = _mm_mul_ps(a.x, b.x);
    result.y = _mm_mul_ps(a.y, b.y);
    return(result);
}

static xvec2_t xvec2_muls(xvec2_t a, __m128 b)
{
    xvec2_t result;
    result.x = _mm_mul_ps(a.x, b);
    result.y = _mm_mul_ps(a.y, b);
    return(result);
}

static xvec2_t xvec2_divs(xvec2_t a, __m128 b)
{
    xvec2_t result;
    result.x = _mm_div_ps(a.x, b);
    result.y = _mm_div_ps(a.y, b);
    return(result);
}

static __m128 xvec2_length_sqr(xvec2_t v)
{
    __m128 result = _mm_add_ps(_mm_mul_ps(v.x, v.x), _mm_mul_ps(v.y, v.y));
    return(result);
}

static __m128 xvec2_length(xvec2_t v)
{
    __m128 result = _mm_sqrt_ps(xvec2_length_sqr(v));
    return(result);
}

static xvec2_t xvec2_normalize(xvec2_t v)
{
    xvec2_t result = xvec2_divs(v, _mm_sqrt_ps(xvec2_length_sqr(v)));
    return(result);
}

static xvec2_t xvec2_lerp(xvec2_t a, xvec2_t b, __m128 t)
{
    xvec2_t result = xvec2_add(xvec2_muls(a, _mm_sub_ps(_mm_set_ps1(1.0f), t)), xvec2_muls(b, t));
    return(result);
}

static xvec2_t xvec2_pow(xvec2_t a, __m128 b)
{
    xvec2_t result;
    result.x = xfast_pow(a.x, b);
    result.y = xfast_pow(a.y, b);
    return(result);
}

static xvec2_t xvec2_min(xvec2_t a, xvec2_t b)
{
    xvec2_t result;
    result.x = _mm_min_ps(a.x, b.x);
    result.y = _mm_min_ps(a.y, b.y);
    return(result);
}

static xvec2_t xvec2_max(xvec2_t a, xvec2_t b)
{
    xvec2_t result;
    result.x = _mm_max_ps(a.x, b.x);
    result.y = _mm_max_ps(a.y, b.y);
    return(result);
}

static xvec3_t xvec3_set(__m128 x, __m128 y, __m128 z)
{
    xvec3_t result = { x, y, z };
    return(result);
}

static xvec3_t xvec3_set1(__m128 v)
{
    xvec3_t result = { v, v, v };
    return(result);
}

static xvec3_t xvec3_setv(vec3_t v)
{
    xvec3_t result = { _mm_set_ps1(v.x), _mm_set_ps1(v.y), _mm_set_ps1(v.z) };
    return(result);
}

static xvec3_t xvec3_add(xvec3_t a, xvec3_t b)
{
    xvec3_t result;
    result.x = _mm_add_ps(a.x, b.x);
    result.y = _mm_add_ps(a.y, b.y);
    result.z = _mm_add_ps(a.z, b.z);
    return(result);
}

static xvec3_t xvec3_sub(xvec3_t a, xvec3_t b)
{
    xvec3_t result;
    result.x = _mm_sub_ps(a.x, b.x);
    result.y = _mm_sub_ps(a.y, b.y);
    result.z = _mm_sub_ps(a.z, b.z);
    return(result);
}

static xvec3_t xvec3_mul(xvec3_t a, xvec3_t b)
{
    xvec3_t result;
    result.x = _mm_mul_ps(a.x, b.x);
    result.y = _mm_mul_ps(a.y, b.y);
    result.z = _mm_mul_ps(a.z, b.z);
    return(result);
}

static xvec3_t xvec3_muls(xvec3_t a, __m128 b)
{
    xvec3_t result;
    result.x = _mm_mul_ps(a.x, b);
    result.y = _mm_mul_ps(a.y, b);
    result.z = _mm_mul_ps(a.z, b);
    return(result);
}

static xvec3_t xvec3_divs(xvec3_t a, __m128 b)
{
    xvec3_t result;
    result.x = _mm_div_ps(a.x, b);
    result.y = _mm_div_ps(a.y, b);
    result.z = _mm_div_ps(a.z, b);
    return(result);
}

static __m128 xvec3_length_sqr(xvec3_t v)
{
    __m128 result = _mm_add_ps(_mm_add_ps(_mm_mul_ps(v.x, v.x), _mm_mul_ps(v.y, v.y)), _mm_mul_ps(v.z, v.z));
    return(result);
}

static __m128 xvec3_length(xvec3_t v)
{
    __m128 result = _mm_sqrt_ps(xvec3_length_sqr(v));
    return(result);
}

static xvec3_t xvec3_normalize(xvec3_t v)
{
    xvec3_t result = xvec3_divs(v, _mm_sqrt_ps(xvec3_length_sqr(v)));
    return(result);
}

static xvec3_t xvec3_lerp(xvec3_t a, xvec3_t b, __m128 t)
{
    xvec3_t result = xvec3_add(xvec3_muls(a, _mm_sub_ps(_mm_set_ps1(1.0f), t)), xvec3_muls(b, t));
    return(result);
}

static xvec3_t xvec3_pow(xvec3_t a, __m128 b)
{
    xvec3_t result;
    result.x = xfast_pow(a.x, b);
    result.y = xfast_pow(a.y, b);
    result.z = xfast_pow(a.z, b);
    return(result);
}

static xvec3_t xvec3_min(xvec3_t a, xvec3_t b)
{
    xvec3_t result;
    result.x = _mm_min_ps(a.x, b.x);
    result.y = _mm_min_ps(a.y, b.y);
    result.z = _mm_min_ps(a.z, b.z);
    return(result);
}

static xvec3_t xvec3_max(xvec3_t a, xvec3_t b)
{
    xvec3_t result;
    result.x = _mm_max_ps(a.x, b.x);
    result.y = _mm_max_ps(a.y, b.y);
    result.z = _mm_max_ps(a.z, b.z);
    return(result);
}

static __m128 xvec3_dot(xvec3_t a, xvec3_t b)
{
    __m128 result = _mm_add_ps(_mm_add_ps(_mm_mul_ps(a.x, b.x), _mm_mul_ps(a.y, b.y)), _mm_mul_ps(a.z, b.z));
    return(result);
}

static xvec3_t xvec3_cross(xvec3_t a, xvec3_t b)
{
    xvec3_t result;
    result.x = _mm_sub_ps(_mm_mul_ps(a.y, b.z), _mm_mul_ps(b.y, a.z));
    result.y = _mm_sub_ps(_mm_mul_ps(a.z, b.x), _mm_mul_ps(b.z, a.x));
    result.z = _mm_sub_ps(_mm_mul_ps(a.x, b.y), _mm_mul_ps(b.x, a.y));
    return(result);
}

static xvec4_t xvec4_set1(__m128 v)
{
    xvec4_t result = { v, v, v, v };
    return(result);
}

static xvec4_t xvec4_set(float x, float y, float z, float w)
{
    xvec4_t result = { _mm_set_ps1(x), _mm_set_ps1(y), _mm_set_ps1(z), _mm_set_ps1(w) };
    return(result);
}

static xvec4_t xvec4_setv(vec4_t v)
{
    xvec4_t result = { _mm_set_ps1(v.x), _mm_set_ps1(v.y), _mm_set_ps1(v.z), _mm_set_ps1(v.w) };
    return(result);
}

static xvec4_t xvec4_add(xvec4_t a, xvec4_t b)
{
    xvec4_t result;
    result.x = _mm_add_ps(a.x, b.x);
    result.y = _mm_add_ps(a.y, b.y);
    result.z = _mm_add_ps(a.z, b.z);
    result.w = _mm_add_ps(a.w, b.w);
    return(result);
}

static xvec4_t xvec4_sub(xvec4_t a, xvec4_t b)
{
    xvec4_t result;
    result.x = _mm_sub_ps(a.x, b.x);
    result.y = _mm_sub_ps(a.y, b.y);
    result.z = _mm_sub_ps(a.z, b.z);
    result.w = _mm_sub_ps(a.w, b.w);
    return(result);
}

static xvec4_t xvec4_mul(xvec4_t a, xvec4_t b)
{
    xvec4_t result;
    result.x = _mm_mul_ps(a.x, b.x);
    result.y = _mm_mul_ps(a.y, b.y);
    result.z = _mm_mul_ps(a.z, b.z);
    result.w = _mm_mul_ps(a.w, b.w);
    return(result);
}

static xvec4_t xvec4_muls(xvec4_t a, __m128 b)
{
    xvec4_t result;
    result.y = _mm_mul_ps(a.y, b);
    result.z = _mm_mul_ps(a.z, b);
    result.x = _mm_mul_ps(a.x, b);
    result.w = _mm_mul_ps(a.w, b);
    return(result);
}

static xvec4_t xvec4_divs(xvec4_t a, __m128 b)
{
    xvec4_t result;
    result.y = _mm_div_ps(a.y, b);
    result.z = _mm_div_ps(a.z, b);
    result.x = _mm_div_ps(a.x, b);
    result.w = _mm_div_ps(a.w, b);
    return(result);
}

static __m128 xvec4_length_sqr(xvec4_t v)
{
    __m128 result = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(v.x, v.x), _mm_mul_ps(v.y, v.y)), _mm_mul_ps(v.z, v.z)), _mm_mul_ps(v.w, v.w));
    return(result);
}
    
static __m128 xvec4_length(xvec4_t v)
{
    __m128 result = _mm_sqrt_ps(xvec4_length_sqr(v));
    return(result);
}

static xvec4_t xvec4_normalize(xvec4_t v)
{
    xvec4_t result = xvec4_divs(v, _mm_sqrt_ps(xvec4_length_sqr(v)));
    return(result);
}

static xvec4_t xvec4_lerp(xvec4_t a, xvec4_t b, __m128 t)
{
    xvec4_t result = xvec4_add(xvec4_muls(a, _mm_sub_ps(_mm_set_ps1(1.0f), t)), xvec4_muls(b, t));
    return(result);
}

static xvec4_t xvec4_pow(xvec4_t a, __m128 b)
{
    xvec4_t result;
    result.x = xfast_pow(a.x, b);
    result.y = xfast_pow(a.y, b);
    result.z = xfast_pow(a.z, b);
    result.w = xfast_pow(a.w, b);
    return(result);
}

static __m128 xvec4_dot(xvec4_t a, xvec4_t b)
{
    __m128 result = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(a.x, b.x), _mm_mul_ps(a.y, b.y)), _mm_mul_ps(a.z, b.z)), _mm_mul_ps(a.w, b.w));
    return(result);
}

static xvec4_t xvec4_cross(xvec4_t a, xvec4_t b)
{
    xvec4_t result;
    result.x = _mm_sub_ps(_mm_mul_ps(a.y, b.z), _mm_mul_ps(b.y, a.z));
    result.y = _mm_sub_ps(_mm_mul_ps(a.z, b.x), _mm_mul_ps(b.z, a.x));
    result.z = _mm_sub_ps(_mm_mul_ps(a.x, b.y), _mm_mul_ps(b.x, a.y));
    result.w = _mm_set_ps1(0.0f);
    return(result);
}

static xvec4_t xvec4_min(xvec4_t a, xvec4_t b)
{
    xvec4_t result;
    result.x = _mm_min_ps(a.x, b.x);
    result.y = _mm_min_ps(a.y, b.y);
    result.z = _mm_min_ps(a.z, b.z);
    result.w = _mm_min_ps(a.w, b.w);
    return(result);
}

static xvec4_t xvec4_max(xvec4_t a, xvec4_t b)
{
    xvec4_t result;
    result.x = _mm_max_ps(a.x, b.x);
    result.y = _mm_max_ps(a.y, b.y);
    result.z = _mm_max_ps(a.z, b.z);
    result.w = _mm_max_ps(a.w, b.w);
    return(result);
}

static xvec3_t xm3x3_mulv(xm3x3_t a, xvec3_t b)
{
    xvec3_t result;
    __m128 x00 = _mm_mul_ps(b.x, a.e[0][0]);
    __m128 x01 = _mm_mul_ps(b.x, a.e[0][1]);
    __m128 x02 = _mm_mul_ps(b.x, a.e[0][2]);
    __m128 y10 = _mm_mul_ps(b.y, a.e[1][0]);
    __m128 y11 = _mm_mul_ps(b.y, a.e[1][1]);
    __m128 y12 = _mm_mul_ps(b.y, a.e[1][2]);
    __m128 z20 = _mm_mul_ps(b.z, a.e[2][0]);
    __m128 z21 = _mm_mul_ps(b.z, a.e[2][1]);
    __m128 z22 = _mm_mul_ps(b.z, a.e[2][2]);
    result.x = _mm_add_ps(_mm_add_ps(x00, y10), z20);
    result.y = _mm_add_ps(_mm_add_ps(x01, y11), z21);
    result.z = _mm_add_ps(_mm_add_ps(x02, y12), z22);
    return(result);
}

static xm4x4_t xm4x4_setm(m4x4_t m)
{
    xm4x4_t result;
    for(uint32_t i = 0; i < 4; i++)
        for(uint32_t j = 0; j < 4; j++)
            result.e[i][j] = _mm_set_ps1(m.e[i][j]);
    return(result);
}

static xvec4_t xm4x4_mulv(xm4x4_t a, xvec4_t b)
{
    xvec4_t result;
    __m128 x00 = _mm_mul_ps(b.x, a.e[0][0]);
    __m128 x01 = _mm_mul_ps(b.x, a.e[0][1]);
    __m128 x02 = _mm_mul_ps(b.x, a.e[0][2]);
    __m128 x03 = _mm_mul_ps(b.x, a.e[0][3]);
    __m128 y10 = _mm_mul_ps(b.y, a.e[1][0]);
    __m128 y11 = _mm_mul_ps(b.y, a.e[1][1]);
    __m128 y12 = _mm_mul_ps(b.y, a.e[1][2]);
    __m128 y13 = _mm_mul_ps(b.y, a.e[1][3]);
    __m128 z20 = _mm_mul_ps(b.z, a.e[2][0]);
    __m128 z21 = _mm_mul_ps(b.z, a.e[2][1]);
    __m128 z22 = _mm_mul_ps(b.z, a.e[2][2]);
    __m128 z23 = _mm_mul_ps(b.z, a.e[2][3]);
    __m128 w30 = _mm_mul_ps(b.w, a.e[3][0]);
    __m128 w31 = _mm_mul_ps(b.w, a.e[3][1]);
    __m128 w32 = _mm_mul_ps(b.w, a.e[3][2]);
    __m128 w33 = _mm_mul_ps(b.w, a.e[3][3]);
    result.x = _mm_add_ps(_mm_add_ps(_mm_add_ps(x00, y10), z20), w30);
    result.y = _mm_add_ps(_mm_add_ps(_mm_add_ps(x01, y11), z21), w31);
    result.z = _mm_add_ps(_mm_add_ps(_mm_add_ps(x02, y12), z22), w32);
    result.w = _mm_add_ps(_mm_add_ps(_mm_add_ps(x03, y13), z23), w33);
    return(result);
}

#ifdef __cplusplus
}
#endif
