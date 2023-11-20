#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commdlg.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "objzero.h"
#include "objzero.c"

#if defined(_DEBUG)

#define SCREEN_SCALE 2
#define SCREEN_W 400
#define SCREEN_H 300

#else

#define SCREEN_SCALE 1
#define SCREEN_W 1440
#define SCREEN_H 900

#endif

#define OFFSET_OF(type, member)(uintptr_t)(&((type *)0)->member)

#include "no_gl.h"

typedef enum menu_t
{
    MENU_OPEN = 1,
    MENU_QUIT,
} menu_t;

typedef struct camera_t
{
    vec3_t pos;
    vec3_t rot;
    float fov;
    float nz;
    float fz;
} camera_t;

typedef struct material_t
{
    vec4_t diffuse_color;
    float shininess;
    texture_t *diffuse_map;
    texture_t *normal_map;
} material_t;

typedef struct mesh_t
{
    uint32_t first_index;
    uint32_t first_vertex;
    uint32_t num_indices;
    uint32_t num_vertices;
    uint32_t material_index;
} mesh_t;

typedef struct model_t
{
    mesh_t *meshes;
    vertex_t *vertices;
    uint32_t *indices;
    uint32_t num_meshes;
    uint32_t num_indices;
    uint32_t num_vertices;
    m4x4_t transform;
} model_t;

typedef struct scene_t
{
    vec3_t bg_color;
    material_t *materials;
    model_t *models;
    uint32_t num_models;
    uint32_t num_materials;
    camera_t cam;
} scene_t;

scene_t scene = { 0 };

static svec4_t blinn_phong_shader(svec3_t pos, svec3_t normal, svec3_t tangent, svec2_t uv)
{
    svec4_t color = get_uniform_vec4(4);
    if(get_texture(0))
        color = svec4_mul(color, sample_texture(get_texture(0), uv));

    svec3_t n = svec3_normalize(normal);
    if(get_texture(1))
    {
        svec3_t nt = svec3_sub(svec3_muls(sample_texture(get_texture(1), uv).xyz, scalar(2.0f)), svec3_set1(scalar(1.0f)));
        svec3_t t = svec3_normalize(sm4x4_mulv(get_model_matrix(), (svec4_t){ tangent.x, tangent.y, tangent.z, scalar(0.0f) }).xyz);
        t = svec3_normalize(svec3_sub(t, svec3_muls(n, svec3_dot(t, n))));
        svec3_t b = svec3_cross(n, t);
        sm3x3_t tbn = { t.x, t.y, t.z, b.x, b.y, b.z, n.x, n.y, n.z };
        n = svec3_normalize(sm3x3_mulv(tbn, nt));
    }

    svec3_t light_dir = get_uniform_vec3(1);
    scalar_t ambient = scalar(0.2f);
    scalar_t diffuse = smax(svec3_dot(light_dir, n), scalar(0.0f));
    color.rgb = svec3_add(svec3_muls(color.rgb, ambient), svec3_muls(color.rgb, diffuse));

    svec3_t light_color = get_uniform_vec3(2);
    svec3_t view_pos = get_uniform_vec3(0);
    svec3_t view_dir = svec3_normalize(svec3_sub(view_pos, pos));
    svec3_t h = svec3_normalize(svec3_add(view_dir, light_dir));
    scalar_t shininess = get_uniform_scalar(5);
    scalar_t specular = spow(smax(svec3_dot(n, h), scalar(0.0f)), shininess);
    color.rgb = svec3_add(color.rgb, svec3_muls(light_color, specular));

    return(color);
}

void find_directory_in_path(char *path, char *result)
{
    char *separator = 0;
    char *str = result;
    while(*path && !(*path == '/' && !*(path + 1)))
    {
        if(*path == '/')
            separator = str;
        *str++ = *path++;
    }
    if(separator)
        *separator = '\0';
    else
        strcpy(result, ".");
}

texture_t *load_texture(char *path)
{
    texture_t *texture = NULL;
    stbi_set_flip_vertically_on_load(true);
    int32_t w, h, comp;
    uint8_t *bytes = stbi_load(path, &w, &h, &comp, 0);
    if(bytes)
    {
        texture = create_texture(w, h, ADDRESS_MODE_REPEAT);
        for(int32_t i = 0; i < (w*h); i++)
        {
            uint32_t r, g, b, a;
            r = g = b = a = 255;

            if(comp == 1)
            {
                r = bytes[i];
                g = bytes[i];
                b = bytes[i];
            }
            else if(comp == 2)
            {
                r = bytes[2*i];
                g = bytes[2*i + 1];
                b = 0;
            }
            else if(comp == 3)
            {
                r = bytes[3*i];
                g = bytes[3*i + 1];
                b = bytes[3*i + 2];
            }
            else if(comp == 4)
            {
                r = bytes[4*i];
                g = bytes[4*i + 1];
                b = bytes[4*i + 2];
                a = bytes[4*i + 3];
            }

            texture->pixels[i] = b | (g << 8) | (r << 16) | (a << 24);
        }
        stbi_image_free(bytes);
    }
    return(texture);
}

model_t *load_model(char *path)
{
    model_t *model = NULL;

    objz_setIndexFormat(OBJZ_INDEX_FORMAT_U32);
    objz_setVertexFormat(sizeof(vertex_t), OFFSET_OF(vertex_t, pos), OFFSET_OF(vertex_t, uv), OFFSET_OF(vertex_t, normal));
    objzModel *obj_model = objz_load(path);

    if(obj_model)
    {
        uint32_t idx = scene.num_models++;
        scene.models = realloc(scene.models, sizeof(model_t)*scene.num_models);
        model = &scene.models[idx];

        memset(model, 0, sizeof(model_t));
        model->transform = m4x4_identity();        
        model->num_indices = obj_model->numIndices;
        model->num_vertices = obj_model->numVertices;
        model->num_meshes = obj_model->numMeshes;

        uint32_t material_offset = scene.num_materials;
        scene.num_materials += obj_model->numMaterials;
        scene.materials = realloc(scene.materials, scene.num_materials*sizeof(material_t));

        model->meshes = malloc(sizeof(mesh_t)*obj_model->numMeshes);
        model->vertices = malloc(sizeof(vertex_t)*obj_model->numVertices);
        memset(model->vertices, 0, sizeof(vertex_t)*obj_model->numVertices);
        model->indices = malloc(sizeof(uint32_t)*obj_model->numIndices);

        for(uint32_t i = 0; i < obj_model->numVertices; i++)
            model->vertices[i] = ((vertex_t *)obj_model->vertices)[i];
        for(uint32_t i = 0; i < obj_model->numIndices; i++)
            model->indices[i] = ((uint32_t *)obj_model->indices)[i];

        char dir[MAX_PATH];
        find_directory_in_path(path, dir);

        for(uint32_t i = 0; i < obj_model->numMaterials; i++)
        {
            objzMaterial *obj_material = &obj_model->materials[i];
            material_t *material = &scene.materials[material_offset + i];
            memset(material, 0, sizeof(material_t));
            material->diffuse_color = (vec4_t){ obj_material->diffuse[0], obj_material->diffuse[1], obj_material->diffuse[2], 1 };
            material->shininess = obj_material->specularExponent;
            if(strlen(obj_material->diffuseTexture) > 0)
            {
                char abs_path[MAX_PATH];
                sprintf(abs_path, "%s/%s", dir, obj_material->diffuseTexture);
                material->diffuse_map = load_texture(abs_path);
            }
            if(strlen(obj_material->bumpTexture) > 0)
            {
                char abs_path[MAX_PATH];
                sprintf(abs_path, "%s/%s", dir, obj_material->diffuseTexture);
                material->normal_map = load_texture(abs_path);
            }
        }

        for(uint32_t i = 0; i < obj_model->numMeshes; i++)
        {
            objzMesh *obj_mesh = &obj_model->meshes[i];
            mesh_t *mesh = &model->meshes[i];
            mesh->first_index = obj_mesh->firstIndex;
            mesh->num_indices = obj_mesh->numIndices;
            mesh->material_index = 0;
            if(obj_mesh->materialIndex >= 0)
                mesh->material_index = material_offset + obj_mesh->materialIndex;
        }

        objz_destroy(obj_model);

        FILE *file = fopen(path, "r");
        if(file)
        {
            uint32_t num = 0;
            float x, y, z, w;
            char line[256];
            while(fgets(line, sizeof(line), file))
            {
                int32_t result = sscanf(line, "# ext.tangent %f %f %f %f", &x, &y, &z, &w);
                if(result == 4)
                {
                    uint32_t i = num++;
                    model->vertices[i].tangent.x = x;
                    model->vertices[i].tangent.y = y;
                    model->vertices[i].tangent.z = z;
                }
            }
            fclose(file);
        }
    }

    return(model);
}

void scale_and_center_scene(void)
{
    vec3_t min_bounds = vec3_set1(FLT_MAX);
    vec3_t max_bounds = vec3_set1(-FLT_MAX);
    for(uint32_t i = 0; i < scene.num_models; i++)
    {
        model_t *model = &scene.models[i];
        for(uint32_t j = 0; j < model->num_vertices; j++)
        {
            min_bounds = vec3_min(min_bounds, model->vertices[j].pos);
            max_bounds = vec3_max(max_bounds, model->vertices[j].pos);
        }
    }

    vec3_t center = vec3_muls(vec3_add(min_bounds, max_bounds), 0.5f);
    vec3_t extent = vec3_muls(vec3_sub(max_bounds, min_bounds), 0.5f);
    float scale = 10.0f/MAX(MAX(extent.x, extent.y), extent.z);
    m4x4_t transform = m4x4_mul(m4x4_scale(vec3_set1(scale)), m4x4_translate(vec3_muls(center, -1.0f)));

    for(uint32_t i = 0; i < scene.num_models; i++)
    {
        model_t *model = &scene.models[i];
        model->transform = m4x4_mul(transform, model->transform);
    }
}

void reset_scene(void)
{
    for(uint32_t i = 0; i < scene.num_models; i++)
    {
        model_t *model = &scene.models[i];
        free(model->meshes);
        free(model->indices);
        free(model->vertices);
    }
    free(scene.models);
    scene.models = NULL;
    scene.num_models = 0;

    for(uint32_t i = 0; i < scene.num_materials; i++)
    {
        material_t *material = &scene.materials[i];
        if(material->diffuse_map)
            free_texture_memory(material->diffuse_map);
        if(material->normal_map)
            free_texture_memory(material->normal_map);
    }
    free(scene.materials);
    scene.materials = NULL;
    scene.num_materials = 0;

    scene.bg_color = (vec3_t){ 0.15f, 0.15f, 0.15f };
    scene.cam = (camera_t)
        {
            .pos = { 0, 0, 35 },
            .rot = { .y = PI },
            .fov = 0.5f*HALF_PI,
            .nz = 0.3f,
            .fz = 150.0f
        };
    scene.num_materials = 1;
    scene.materials = malloc(sizeof(material_t));
    scene.materials[0] = (material_t)
        {
            .diffuse_color = { 0.5f, 0.5f, 0.5f, 1.0f },
            .shininess = 32.0f,
        };
}

static m4x4_t camera_view(camera_t *cam)
{
    m4x4_t tm = m4x4_rotate(cam->rot);
    vec3_t forward = m4x4_mulv(tm, (vec4_t){ 0, 0, 1, 0 }).xyz;
    vec3_t up = m4x4_mulv(tm, (vec4_t){ 0, 1, 0, 0 }).xyz;
    m4x4_t view = m4x4_look_at(cam->pos, vec3_add(cam->pos, forward), up);
    return(view);
}

static m4x4_t camera_projection(camera_t *cam, vec2_t screen_size)
{
    m4x4_t proj = m4x4_perspective(cam->fov, screen_size.x/screen_size.y, cam->nz, cam->fz);
    return(proj);
}

static void convert_backslashes(char *path)
{
    while(*path)
        *path++ = *path == '\\' ? '/' : *path;
}

static LRESULT CALLBACK window_callback(HWND window, UINT msg, WPARAM wparam, LPARAM lparam)
{
    LRESULT result = 0;
    switch(msg)
    {
        case WM_CLOSE:
        {
            DestroyWindow(window);
            PostQuitMessage(0);
            break;
        }
        case WM_COMMAND:
        {
            switch(LOWORD(wparam))
            {
                case MENU_OPEN:
                {
                    char path[MAX_PATH] = { 0 };
                    OPENFILENAMEA open_file =
                        {
                            .lStructSize = sizeof(OPENFILENAMEA),
                            .hwndOwner = window,
                            .lpstrFile = path,
                            .nMaxFile = sizeof(path),
                            .lpstrFilter = "Wavefront Files (*.obj)\0*.obj\0",
                            .nFilterIndex = 1,
                            .Flags = (OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR),
                        };
                    if(GetOpenFileName(&open_file))
                    {
                        convert_backslashes(path);
                        reset_scene();
                        load_model(path);
                        scale_and_center_scene();
                    }
                    break;
                }
                case MENU_QUIT:
                {
                    SendMessage(window, WM_CLOSE, 0, 0);
                    break;
                }
            }
            break;
        }
        default:
        {
            result = DefWindowProc(window, msg, wparam, lparam);
            break;
        }
    }
    return(result);
}

int main(int argc, char *argsv[])
{
    WNDCLASS window_class =
        {
            .style = (CS_HREDRAW | CS_VREDRAW | CS_OWNDC),
            .lpfnWndProc = window_callback,
            .hInstance = GetModuleHandle(NULL),
            .hIcon = LoadIcon(NULL, IDI_APPLICATION),
            .hCursor = LoadCursor(NULL, IDC_ARROW),
            .hbrBackground = GetStockObject(NULL_BRUSH),
            .lpszClassName = "window_class",
        };

    if(RegisterClass(&window_class))
    {
        LONG window_w = SCREEN_SCALE*SCREEN_W;
        LONG window_h = SCREEN_SCALE*SCREEN_H;
        RECT r = { .right = window_w, .bottom = window_h };
        DWORD window_style = (WS_OVERLAPPEDWINDOW & ~(WS_THICKFRAME | WS_MINIMIZEBOX | WS_MAXIMIZEBOX));
        AdjustWindowRectEx(&r, window_style, true, 0);
        LONG w = (r.right - r.left);
        LONG h = (r.bottom - r.top);
        HWND window = CreateWindowEx(0, window_class.lpszClassName, "noGL Viewer", window_style, CW_USEDEFAULT, CW_USEDEFAULT, w, h, NULL, NULL, window_class.hInstance, NULL);

        if(window)
        {
            HMENU menu = CreateMenu();
            AppendMenu(menu, MF_STRING, MENU_OPEN, "Open...");
            AppendMenu(menu, MF_SEPARATOR, 0, NULL);
            AppendMenu(menu, MF_STRING, MENU_QUIT, "Quit");

            HMENU menubar = CreateMenu();
            AppendMenu(menubar, MF_POPUP, (UINT_PTR)menu, "File");
            SetMenu(window, menubar);

            ShowWindow(window, SW_SHOWNORMAL);
            HDC context = GetDC(window);

            BITMAPINFO bmpi =
                {
                    .bmiHeader.biSize = sizeof(BITMAPINFOHEADER),
                    .bmiHeader.biWidth = SCREEN_W,
                    .bmiHeader.biHeight = -SCREEN_H,
                    .bmiHeader.biPlanes = 1,
                    .bmiHeader.biBitCount = 32,
                    .bmiHeader.biCompression = BI_RGB,
                };

            reset_scene();

            setup_gl();

            framebuffer_t *framebuffer = create_framebuffer(SCREEN_W, SCREEN_H);
            set_framebuffer(framebuffer);

            LARGE_INTEGER frequency;
            QueryPerformanceFrequency(&frequency);
            LARGE_INTEGER ticks;
            QueryPerformanceCounter(&ticks);

            vec2_t mouse_pos = { 0, 0 };

            bool running = true;
            while(running)
            {
                MSG msg;
                while(PeekMessage(&msg, 0, 0, 0, PM_REMOVE))
                {
                    running = (msg.message != WM_QUIT);
                    if(!running)
                        break;
                    TranslateMessage(&msg);
                    DispatchMessage(&msg);
                }

                LARGE_INTEGER old_ticks = ticks;
                QueryPerformanceCounter(&ticks);
                double elapsed = ((double)(ticks.QuadPart - old_ticks.QuadPart)/(double)frequency.QuadPart);

                char buffer[128];
                sprintf(buffer, "noGL Viewer - %.2lfms", (elapsed*1000.0));
                SetWindowTextA(window, buffer);

                POINT cursor_pos;
                GetCursorPos(&cursor_pos);
                ScreenToClient(window, &cursor_pos);

                vec2_t old_mouse_pos = mouse_pos;
                mouse_pos = (vec2_t){ (float)cursor_pos.x, (float)cursor_pos.y };
                vec2_t mouse_delta = vec2_sub(mouse_pos, old_mouse_pos);

                camera_t *cam = &scene.cam;

                vec3_t move = { 0 };
                m4x4_t tm = m4x4_rotate(cam->rot);
                vec3_t forward = m4x4_mulv(tm, (vec4_t){ 0, 0, 1, 0 }).xyz;
                vec3_t right = m4x4_mulv(tm, (vec4_t){ 1, 0, 0, 0 }).xyz;

                if(GetActiveWindow() == window)
                {
                    if(((GetKeyState(VK_RBUTTON) & 0x8000) != 0))
                    {
                        cam->rot.x += mouse_delta.y*0.003f;
                        cam->rot.y += mouse_delta.x*0.003f;
                        cam->rot.x = clamp(cam->rot.x, -HALF_PI, HALF_PI);
                    }

                    if(((GetKeyState('D') & 0x8000) != 0))
                        move = vec3_add(move, right);
                    if(((GetKeyState('A') & 0x8000) != 0))
                        move = vec3_sub(move, right);
                    if(((GetKeyState('W') & 0x8000) != 0))
                        move = vec3_add(move, forward);
                    if(((GetKeyState('S') & 0x8000) != 0))
                        move = vec3_sub(move, forward);
                }

                if(vec3_length_sqr(move) > 0)
                {
                    float speed = (float)elapsed*15.0f;
                    cam->pos = vec3_add(cam->pos, vec3_muls(vec3_normalize(move), speed));
                }

                vec3_t light_color = { 0.3f, 0.2f, 0.15f };
                vec3_t light_dir = vec3_muls(forward, -1.f);

                vec2_t screen_size = { (float)SCREEN_W, (float)SCREEN_H };

                set_view_matrix(camera_view(cam));
                set_projection_matrix(camera_projection(cam, screen_size));

                clear_color(framebuffer, (vec4_t){ scene.bg_color.r, scene.bg_color.r, scene.bg_color.r, 0.0f });
                clear_depth(framebuffer, 1.0f);

                set_fragment_shader(&blinn_phong_shader);
                set_uniform_vec3(0, cam->pos);
                set_uniform_vec3(1, light_dir);
                set_uniform_vec3(2, light_color);

                for(uint32_t i = 0; i < scene.num_models; i++)
                {
                    model_t *model = &scene.models[i];
                    set_model_matrix(model->transform);
                    for(uint32_t j = 0; j < model->num_meshes; j++)
                    {
                        mesh_t *mesh = &model->meshes[j];
                        material_t *material = &scene.materials[mesh->material_index];
                        set_uniform_vec4(4, material->diffuse_color);
                        set_uniform_scalar(5, material->shininess);
                        set_texture(0, material->diffuse_map);
                        set_texture(1, material->normal_map);
                        draw(model->vertices, model->indices + mesh->first_index, mesh->num_indices);
                    }
                }

                StretchDIBits(context, 0, 0, window_w, window_h, 0, 0, SCREEN_W, SCREEN_H, framebuffer->color, &bmpi, DIB_RGB_COLORS, SRCCOPY);
            }
        }
    }
    return(0);
}
