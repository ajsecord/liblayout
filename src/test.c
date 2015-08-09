
#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>   /* For timing with gettimeofday() */
#include <string.h>

#include <layout/layout.h>

#include "random/random.h"

typedef struct {
    lay_coord_t x;          /**< The x-coordinate of the "left" edge (or, smallest x value). */
    lay_coord_t y;          /**< The y-coordinate of the "top" edge (or, smallest y value). */
    lay_extent_t width;     /**< The width. */
    lay_extent_t height;    /**< The height. */
    int fixed;              /**< The fixed flag. */
} rect;

typedef struct {
    int size;               /**< The number of rectangles in this array. */
    rect* items;            /**< The rectangles. */
} rect_array;

/* Distribution of rectangles */
static int default_num_rects = 200;
static lay_real_t default_extent_mean = 40;
static lay_real_t default_extent_stddev = 0;
static lay_real_t default_extent_min = 5;

/* Energy weights */
static const lay_real_t overlap_weight = 5e2;
static const lay_real_t edge_weight = 1e2;
static const lay_real_t center_weight = 0;
static const lay_real_t original_pos_weight = 1;

#if 0
/* Simulated annealing */
static const lay_real_t max_temp = 10;
static const lay_real_t temp_rate = 1e-4;
static lay_real_t cur_temp = -1;
#endif

/* Drawing */
static const GLfloat rect_color[4]       = { 0.2, 0.6, 0.7, 0.75 };
static const GLfloat overlap_color[4]    = { 0.8, 0.6, 0.7, 0.75 };
static const GLfloat fixed_color[4]      = { 0.5, 0.5, 0.5, 0.75 };
static const GLfloat line_color[4]       = { 0.2, 0.6, 0.9, 0.75 };
static const GLfloat grad_color[4]       = { 0.9, 0.2, 0.2, 1    };
static const GLfloat hover_edge_color[4] = {   0,   0,   0, 0.5  };
static const GLfloat drag_edge_color[4]  = {   0,   0,   0, 1    };

/* Application data */
static const int num_steps_per_idle = 1;
static rect_array *rects = NULL;
static long seed = 0;
static int screen_width = -1, screen_height = -1;
static int animating = 0;

static int sel_rect = -1;
static int drag = 0, drag_start_x, drag_start_y;

typedef enum { SIM_ANNEAL, MACOPT, NUM_METHODS } method_type;
static method_type method = MACOPT;

/* Function declarations */
static void idle(void);

static rect_array* rect_array_create(const int size) {
    rect_array* result = malloc(sizeof(rect_array));
    assert(result);
    
    result->size = size;
    result->items = malloc(sizeof(rect) * size);
    
    return result;
}

static void rect_array_copy(const rect_array const* in, rect_array* out) {
    assert(in && out);
    assert(in->size == out->size);
    memcpy(out->items, in->items, in->size * sizeof(rect));
}

static void rect_array_destroy(rect_array* array) {
    if (!array)
        return;
    
    assert((array->size == 0) == (array->items == NULL));
    assert(array->size >= 0);
    
    if (array->size > 0) {
        assert(array->items);
        free(array->items);
    }
    free(array);
}

int point_in_rect(const lay_coord_t x, const lay_coord_t y, const rect const* r) {
    assert(r);
    return (x >= r->x && x <= (r->x + r->width) && 
            y >= r->y && y <= (r->y + r->height));
}


static int calc_num_dof(const rect_array const* r) {
    int i, count = 0;
    assert(r && r->size >= 0);
    for (i = 0; i < r->size; ++i) {
        count += (r->items[i].fixed == 0);
    }
    
    return 2 * count;
}

static void copy_rect_pos_to_array(const rect_array const* r, lay_real_t* array) {
    int i, j;
    assert(r && array);
    for (i = 0, j = 0; i < r->size; ++i) {
        if (!r->items[i].fixed) {
            array[2*j]   = r->items[i].x;
            array[2*j+1] = r->items[i].y;
            ++j;
        }
    }
    assert(j >= 0 && j <= r->size);
}

static void copy_array_to_rect_pos(const lay_real_t const* array, rect_array* r) {
    int i, j;
    assert(r && array);
    for (i = 0, j = 0; i < r->size; ++i) {
        if (!r->items[i].fixed) {
            r->items[i].x = array[2*j];
            r->items[i].y = array[2*j+1];
            ++j;
        }
    }
    assert(j >= 0 && j <= r->size);
}

#if 0
/* Simulated annealing functions */
static void sa_neighbor(void* context, const void const* input, const lay_real_t cur_temp, void* output) {
    long* rand_seed;
    rect_array in, out;
    int i;
    
    assert(context && input && output);
    rand_seed = (long*) context;
    in = *(const rect_array const*) input;
    out = *(rect_array*) output;
    assert(in.size == out.size);
    
    for (i = 0; i < in.size; ++i) {
        out.items[i].x = in.items[i].x + (uniform_dev(rand_seed) < 0.5 ? -1 : 1) * (lay_coord_t) (2 * uniform_dev(rand_seed));
        out.items[i].y = in.items[i].y + (uniform_dev(rand_seed) < 0.5 ? -1 : 1) * (lay_coord_t) (2 * uniform_dev(rand_seed));
        out.items[i].width = in.items[i].width;
        out.items[i].height = in.items[i].height;
    }
}

static lay_real_t sa_eval(void* context, const void const* input) {
    return eval(input, NULL);
}

static lay_real_t sa_reduce_temp(void* context, const lay_real_t cur_temp) {
    return (1 - temp_rate) * cur_temp;
}

/* Initialization functions */
static void init_point() {
    point_ops.create = create_point;
    point_ops.copy = copy_point;
    point_ops.num_dof = point_num_dof;
    point_ops.destroy = destroy_point;
}

static void init_sa() {
    sa_ops.neighbor = sa_neighbor;
    sa_ops.eval = sa_eval;
    sa_ops.reduce_temp = sa_reduce_temp;
    
    cur_temp = max_temp;
}

#endif

static void init_rects(const int num, const lay_real_t mean, const lay_real_t stddev, const lay_real_t min) {
    int i;
    lay_real_t l;
    
    rect_array_destroy(rects);
    rects = rect_array_create(num);
    
    for (i = 0; i < num; ++i) {
        rects->items[i].x = screen_width * rng_uniform_dev(&seed);
        rects->items[i].y = screen_height * rng_uniform_dev(&seed);
        
        l = mean + stddev * rng_gauss_dev(&seed);
        rects->items[i].width = (l >= min ? l : min);

        l = mean + stddev * rng_gauss_dev(&seed);
        rects->items[i].height = (l >= min ? l : min);
   
        rects->items[i].fixed = 0;
    }
}


static void speed_test(const int trials) {
#if 0
    struct timeval start, stop;
    int i, j;
    lay_real_t t;
    
    gettimeofday(&start, NULL);
    for (i = 0; i < trials; ++i) {
        for (j = 0; j < rects->size; ++j)
            lay_any_overlap(j, *rects);
    }
    gettimeofday(&stop, NULL);
    
    t = (stop.tv_sec - start.tv_sec) + 1e-6 * (stop.tv_usec - start.tv_usec);
    printf("Per-rectangle time (%i trials, %i rectangles): %g s\n", 
           trials, rects->size, t / trials);
#endif
}

static void init_opengl(void) {
    glClearColor(1, 1, 1, 0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

static void animate(const int go_daddy) {
    animating = go_daddy;
    glutIdleFunc(animating ? idle : NULL);
}

#if 0
static int step_sa(const int steps) {
    lay_real_t energy;
    const int iters = lay_sim_anneal(&seed, sa_ops, point_ops, rects, &cur_temp, 0.0, 100.0, steps, &seed);
    if (iters != steps) 
        fprintf(stderr, "Warning, number of iterations mismatch, requested %i, got %i\n", steps, iters);
   
    energy = sa_eval(NULL, rects);
    printf("temp: %g, energy: %g\n", cur_temp, energy);
    
    return (energy == 0);
}
#endif

static int step_mo(const int steps) {
    lay_statep state = NULL;
    
    state = lay_create_state();
    
    lay_register_rects(state, 
                       &(rects->items[0].x),     sizeof(rect), 
                       &(rects->items[0].width), sizeof(rect),
                       rects->size);
                 
    lay_set_overlap_weight(state, overlap_weight);
    lay_set_center_weight(state, center_weight);
    lay_set_edge_weight(state, edge_weight);
    lay_set_orig_pos_weight(state, original_pos_weight);
    
    lay_optimize(state);
    
    lay_destroy_state(state);

    return 1;
}

static void step(const int steps) {
    int done = 0;
    assert(method >= 0 && method < NUM_METHODS);
    switch(method) {
        /*case SIM_ANNEAL:   done = step_sa(steps); break;*/
        case MACOPT:       done = step_mo(steps); break;
        default:
            assert(!"Invalid method");
            break;
    }
    
    if (done) 
        animate(0);
}

static void draw_rects(const rect_array rects) {
    int i;
    const rect *r;
    
    for (i = 0; i < rects.size; ++i) {
        
        if (rects.items[i].fixed) 
            glColor4fv(fixed_color);
        //else if (lay_any_overlap(i, rects))
        //    glColor4fv(overlap_color);
        else 
            glColor4fv(rect_color);
        
        r = rects.items + i;
        glRectf(r->x, r->y, r->x + r->width, r->y + r->height);

        if (sel_rect == i) {
            if (drag) 
                glColor4fv(drag_edge_color);
            else
                glColor4fv(hover_edge_color);
            glLineWidth(1);
            glBegin(GL_LINE_LOOP);
                glVertex2f(r->x,            r->y);
                glVertex2f(r->x + r->width, r->y);
                glVertex2f(r->x + r->width, r->y + r->height);
                glVertex2f(r->x,            r->y + r->height);
            glEnd();
        }
    }
}

static int find_rect(const int x, const int y) {
    int i, index = -1;
    for (i = rects->size - 1; i >= 0; --i) {
        if (point_in_rect(x, screen_height - y, rects->items + i)) {
            index = i;
            break;
        }
    }
    
    assert(index == -1 || (index >= 0 && index < rects->size));
    return index;
}

static void mouse(int button, int state, int x, int y) {
    int modifiers;
    
    sel_rect = find_rect(x,y);
    if (sel_rect == -1)
        return;
    assert(sel_rect >= 0 && sel_rect < rects->size);
    
    modifiers = glutGetModifiers();
    
    if (button == GLUT_LEFT_BUTTON) {
        /* Plain left-click */
        if (!modifiers) {
            if (state == GLUT_DOWN) {
                drag = 1;
                rects->items[sel_rect].fixed = 1;
                drag_start_x = x;
                drag_start_y = y;
                
            } else {
                drag = 0;
                rects->items[sel_rect].fixed = 0;
            }
        
        /* Shift left-click */
        } else if ((modifiers & GLUT_ACTIVE_SHIFT) && (state == GLUT_DOWN)) {
            rects->items[sel_rect].fixed = !rects->items[sel_rect].fixed;
        }
    }
    
    glutPostRedisplay();
}

static void passive_motion(int x, int y) {
    const int index = find_rect(x, y);
    
    if (sel_rect != index) {
        sel_rect = index;
        glutPostRedisplay();
    }
}

static void motion(int x, int y) {
    if (!drag)
        return;
    
    assert(sel_rect >= 0 && sel_rect < rects->size);
    
    rects->items[sel_rect].x += x - drag_start_x;
    rects->items[sel_rect].y += (screen_height - y) - (screen_height - drag_start_y);
    
    drag_start_x = x;
    drag_start_y = y;
    
    animate(1);
    
    glutPostRedisplay();
}

 static void idle(void) {
    step(num_steps_per_idle);
    glutPostRedisplay();
}

static void display(void) {
    if (rects == NULL) {
#if 1
        init_rects(default_num_rects, 
                   default_extent_mean, default_extent_stddev, default_extent_min);
#else 
        init_rects(2, default_extent_mean, default_extent_stddev, default_extent_min);
        rects->items[0].x = screen_width/2.0 + 1;
        rects->items[0].y = screen_width/2.0;
#endif
    }
    
    glClear(GL_COLOR_BUFFER_BIT);
    
    draw_rects(*rects);
    
    glutSwapBuffers();
    
    glutReportErrors();
}

static void reshape(int width, int height) {
    screen_width = width;
    screen_height = height;
    
    glViewport(0, 0, width, height);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, screen_width, 0, screen_height);
    glMatrixMode(GL_MODELVIEW);
}

static void menu(int value) {
    switch (value) {
        case 'q':
        case 'Q':
            exit(0);
            break;
            
        case 'i':
            init_rects(rects->size,
                       default_extent_mean, default_extent_stddev, default_extent_min);
            /*cur_temp = max_temp;*/
            glutPostRedisplay();
            break;
            
        case ' ':
            animate(!animating);
        break;
            
        case '+':
            default_num_rects *= 2;
            init_rects(default_num_rects, 
                       default_extent_mean, default_extent_stddev, default_extent_min);
            glutPostRedisplay();
            break;
            
        case '-':
            default_num_rects = (default_num_rects/2 > 1 ? default_num_rects : 1);
            init_rects(default_num_rects, 
                       default_extent_mean, default_extent_stddev, default_extent_min);
            glutPostRedisplay();
            break;
            
        case 's':
            step(500);
            glutPostRedisplay();
            break;
            
        case 'm':
            method = (method + 1) % NUM_METHODS;
            switch (method) {
                case SIM_ANNEAL:   printf("Simulated annealing\n"); break;
                case MACOPT:       printf("Conjugated gradient\n"); break;
                default:           assert(!"Invalid method"); break;
            }
            break;
            
        case 't':
            speed_test(100);
            break;
            
        default:
            break;
    }
}

static void keyboard(unsigned char c, int x, int y) {
    menu(c);
}

static void init_menu() {
    glutCreateMenu(menu);
    glutAddMenuEntry("Quit 'q'", (int)'q');
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

int main(int argc, char* argv[]) {
    seed = rng_gen_random_seed();
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
  
    glutInitWindowSize(500, 500);
    glutCreateWindow("Layout optimization");

    printf("Open GL vendor:   %s\n", glGetString(GL_VENDOR));
    printf("Open GL renderer: %s\n", glGetString(GL_RENDERER));
    printf("Open GL version:  %s\n", glGetString(GL_VERSION));
   
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutPassiveMotionFunc(passive_motion);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    
    init_menu();
    init_opengl();
#if 0
    init_point();
    init_sa();
#endif
    glutMainLoop();
    return 0;
}
