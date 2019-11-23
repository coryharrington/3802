#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "matrix.h"

#define     NUM_THREADS     2

typedef struct {
    unsigned int id;
    TMatrix *m, *n, *t;
} thread_arg_t;

/* the main function of threads */
static void * thread_main(void * p_arg)
{
    thread_arg_t *arg = p_arg;
    for (unsigned int i = arg->id; i < arg->m->nrows; i += NUM_THREADS)
        for (unsigned int j = 0; j < arg->m->ncols; j++)
            arg->t->data[i][j] = arg->m->data[i][j] + arg->n->data[i][j];
    return NULL;
}

/* Return the sum of two matrices. The result is in a newly creaed matrix. 
 *
 * If a pthread function fails, report error and exit. 
 * Return NULL if something else is wrong.
 *
 * Similar to addMatrix, but this function uses 2 threads.
 */
TMatrix * addMatrix_thread(TMatrix *m, TMatrix *n)
{
    if (    m == NULL || n == NULL
         || m->nrows != n->nrows || m->ncols != n->ncols )
        return NULL;

    TMatrix * t = newMatrix(m->nrows, m->ncols);
    if (t == NULL)
        return t;

    pthread_t thread_ids[NUM_THREADS];
    thread_arg_t thread_args[NUM_THREADS];
    // Prepare arguments, create threads, check error
    for(int i = 0; i < NUM_THREADS; i++)
    {
        thread_args[i].id = i;
        thread_args[i].m = m;
        thread_args[i].n = n;
        thread_args[i].t = t;
        int rv = pthread_create(&thread_ids[i], NULL, thread_main, &thread_args[i]);
        if(rv)
        {
            fprintf(stderr, "ERROR: pthread_create() returned %d\n", rv);
            exit(1);
        }
    }
    // Join threads, check error
    for(int i = 0; i < NUM_THREADS; i++)
    {
        int rv = pthread_join(thread_ids[i], NULL);
        if(rv)
        {
            fprintf(stderr, "ERROR: pthread_join() returned %d\n", rv);
            exit(1);
        }
    }

    return t;
}
