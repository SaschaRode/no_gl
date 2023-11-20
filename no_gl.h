#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#if !defined(_DEBUG)
#define FAST_MODE // Enable fast mode to use SIMD instructions in shaders
#endif

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#define CONCAT_IMPL(a, b) a##b
#define CONCAT(a, b) CONCAT_IMPL(a, b)
#define PAD(n) char CONCAT(_padding_, __LINE__)[n]

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define CLAMP(v, a, b) MIN(MAX((v), (a)), (b))

#if defined(_DEBUG)
#include <assert.h>
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)
#endif

#if defined(_WIN32)

#include <intrin.h>

#define ALIGNED(bytes) __declspec(align(bytes))

#if defined(EXPORT)
#define API __declspec(dllexport)
#else
#define API __declspec(dllimport)
#endif

#else

#define ALIGNED(bytes)
#define API

#endif

#include "no_gl_math.h"

typedef svec4_t fragment_shader_t(svec3_t pos, svec3_t normal, svec3_t tangent, svec2_t uv);

typedef enum address_mode_t
{
    ADDRESS_MODE_REPEAT,
    ADDRESS_MODE_CLAMP,
} address_mode_t;

typedef struct vertex_t
{
    vec3_t pos;
    vec3_t normal;
    vec3_t tangent;
    vec2_t uv;
} vertex_t;

typedef struct texture_t
{
    int32_t w, h;
    uint32_t *pixels;
    address_mode_t address_mode;
} texture_t;

typedef struct framebuffer_t
{
    int32_t w, h;
    uint32_t *color;
    float *depth;
} framebuffer_t;

// initialization and cleanup

API void setup_gl(void);

API framebuffer_t *create_framebuffer(int32_t w, int32_t h);
API texture_t *create_texture(int32_t w, int32_t h, address_mode_t address_mode);

API void free_framebuffer_memory(framebuffer_t *framebuffer);
API void free_texture_memory(texture_t *texture);

// pipeline setup

API void clear_depth(framebuffer_t* framebuffer, float depth);
API void clear_color(framebuffer_t *framebuffer, vec4_t color);

API void set_framebuffer(framebuffer_t *framebuffer);
API void set_texture(uint32_t idx, texture_t *texture);
API void set_fragment_shader(fragment_shader_t *shader);

API void set_model_matrix(m4x4_t matrix);
API void set_view_matrix(m4x4_t matrix);
API void set_projection_matrix(m4x4_t matrix);

API void set_uniform_scalar(uint32_t idx, float value);
API void set_uniform_vec2(uint32_t idx, vec2_t value);
API void set_uniform_vec3(uint32_t idx, vec3_t value);
API void set_uniform_vec4(uint32_t idx, vec4_t value);

// drawing

API void draw(vertex_t *vertices, uint32_t *indices, uint32_t num_indices);

// shader accessors

API texture_t *get_texture(uint32_t idx);
API svec4_t sample_texture(texture_t *texture, svec2_t uv);

API sm4x4_t get_model_matrix(void);
API sm4x4_t get_view_matrix(void);
API sm4x4_t get_projection_matrix(void);

API scalar_t get_uniform_scalar(uint32_t idx);
API svec2_t get_uniform_vec2(uint32_t idx);
API svec3_t get_uniform_vec3(uint32_t idx);
API svec4_t get_uniform_vec4(uint32_t idx);

#ifdef __cplusplus
}
#endif
