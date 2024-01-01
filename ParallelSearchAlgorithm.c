#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct synchronize_structure {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int counter;
    int limit;
} synchronize_type;

synchronize_type *synchronize_all_vars;

void synchronize(synchronize_type *synchronize_vars) {
    pthread_mutex_lock(&synchronize_vars->mutex);
    synchronize_vars->counter++;
    if (synchronize_vars->counter == synchronize_vars->limit) {
        synchronize_vars->counter = 0;
        pthread_cond_broadcast(&synchronize_vars->cond);
    }
    else {
        pthread_cond_wait(&synchronize_vars->cond, &synchronize_vars->mutex);
    }
    pthread_mutex_unlock(&synchronize_vars->mutex);
}

// Radna funkcija dretvene funkcije.
typedef void (*working_function_t)(int, void*);
working_function_t function;
// Argument radne funkcije dretvene funkcije.
static void* argument;
// Argument dretvene funkcije.
struct thread_data {
    int rank;
};
struct thread_data **wd;
int stopped;

// Dretvena funkcija koja ceka.
static void *waiting_thread_function(void *arg) {
    struct thread_data *thread_data = (struct thread_data*) arg;
    const int rank = thread_data->rank;
    for (;;) {
        synchronize(synchronize_all_vars);
        if (stopped) {
            break;
        }
        else {
            function(rank, argument);
            synchronize(synchronize_all_vars);
        }
    }
    return NULL;
}

// Argument za zadatak 1.
typedef struct parallel_search_structure {
    int num_threads;
    int n;
    double y;
    double *x;
    int *lim;
    int *c;
    int *q;
    int *found;
    int *res_i;
} parallel_search_type;
parallel_search_type *parallel_search_arg;

// Radna funkcija
void working_function_for_task_1(int rank, void *arg) {

    parallel_search_type *task_args = (parallel_search_type *) arg;

    if (rank == 0) {
        task_args->lim[0] = 0;
        task_args->lim[1] = task_args->n + 1;
        task_args->c[0] = 0;
        task_args->c[task_args->num_threads + 1] = 1;
    }

    synchronize(synchronize_all_vars);

    while ((task_args->lim[1] - task_args->lim[0]) > task_args->num_threads && task_args->found[0] == 0) {

        if (rank == 0) {
            task_args->q[0] = task_args->lim[0];
            task_args->q[task_args->num_threads + 1] = task_args->lim[1];
        } 

        synchronize(synchronize_all_vars);

        task_args->q[rank] = task_args->lim[0] + rank * (task_args->lim[1] - task_args->lim[0]) / (task_args->num_threads + 1);

        double k;
        if (task_args->q[rank] == 0) {
            k = -1e308;
        } else if (task_args->q[rank] >= task_args->n + 1) {
            k = 1e308;
        } else {
            k = task_args->x[task_args->q[rank] - 1];
        }

        if (task_args->y == k) {
            int i = task_args->q[rank] - 1;
            task_args->found[0] = 1;
            task_args->res_i[0] = i;
        } 
        else {
            if (task_args->found[0] == 0) {
                if (task_args->y > k) {
                    task_args->c[rank] = 0;
                }
                else {
                    task_args->c[rank] = 1;
                }
            }
        }

        synchronize(synchronize_all_vars);

        if (task_args->c[rank] < task_args->c[rank + 1] && task_args->found[0] == 0) {
            task_args->lim[0] = task_args->q[rank];
            task_args->lim[1] = task_args->q[rank + 1];
        }

        synchronize(synchronize_all_vars);

        if (rank == 0 && task_args->c[0] < task_args->c[1] && task_args->found[0] == 0) {
            task_args->lim[0] = task_args->q[0];
            task_args->lim[1] = task_args->q[1];
        }

        synchronize(synchronize_all_vars);
    }

    synchronize(synchronize_all_vars);

    if (rank <= task_args->lim[1] - task_args->lim[0]) {

        double k;
        if (task_args->lim[0] + rank == 0) {
            k = -1e308;
        } else if (task_args->lim[0] + rank >= task_args->n + 1) {
            k = 1e308;
        } else {
            k = task_args->x[task_args->lim[0] + rank - 1];
        }
        
        synchronize(synchronize_all_vars);

        if (task_args->y == k) {
            int i = task_args->lim[0] + rank - 1;
            task_args->found[0] = 1;
            task_args->res_i[0] = i;
        }
        
        synchronize(synchronize_all_vars);

        if (task_args->found[0] == 0) {
            if (task_args->y > k) {
                task_args->c[rank] = 0;
            } 
            else {
                task_args->c[rank] = 1;
            }
        }  

        synchronize(synchronize_all_vars);

        if (task_args->c[rank - 1] < task_args->c[rank] && task_args->found[0] == 0) {
            int i = task_args->lim[0] + rank - 2;
            task_args->found[0] = 1;
            task_args->res_i[0] = i;
        }

        synchronize(synchronize_all_vars);
    }

}

int main() {
    int p = 3;
    int max_num_threads = p + 1;

    synchronize_all_vars = (synchronize_type *)malloc(sizeof(synchronize_type));
    pthread_mutex_init(&synchronize_all_vars->mutex, 0);
    pthread_cond_init(&synchronize_all_vars->cond, 0);
    synchronize_all_vars->limit = max_num_threads;
    synchronize_all_vars->counter = 0;

    parallel_search_arg = (parallel_search_type*) malloc(sizeof(parallel_search_type));
    
    pthread_t *threads = (pthread_t*) malloc(p * sizeof(pthread_t));
    wd = (struct thread_data**) malloc(p * sizeof(struct thread_data*));

    for (int i = 0; i < p; i++ ) {
        wd[i] = (struct thread_data*) malloc(sizeof(struct thread_data));
    }

    stopped = 0;
    for (int i = 0; i < p; i++ ) {
        wd[i]->rank = i + 1; // ako glavna dretva isto izvodi radnu funkciju
        pthread_create(&threads[i], 0, waiting_thread_function, (void*) wd[i]);
    }

    double values[] = {50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0};  
    int size = sizeof(values) / sizeof(values[0]);

    // Izvrednjavanje strukture argumenta radne funkcije za zadatak 1.
    parallel_search_arg->num_threads = p;
    parallel_search_arg->n = size;
    parallel_search_arg->y = 370.0;

    parallel_search_arg->x = (double *)malloc(size * sizeof(double)); 
    
    for (int i = 0; i < size; i++) {
        printf("x[%d] = %f\n", i, values[i]);
    }

    printf("\ny = %f\n\n", parallel_search_arg->y);

    memcpy(parallel_search_arg->x, values, sizeof(values));

    parallel_search_arg->lim = (int *)malloc(2 * sizeof(int));
    parallel_search_arg->c = (int *)malloc((parallel_search_arg->n + 1) * sizeof(int));
    parallel_search_arg->q = (int *)malloc((parallel_search_arg->n + 1) * sizeof(int));
    parallel_search_arg->found = (int *)malloc(1 * sizeof(int));
    parallel_search_arg->res_i = (int *)malloc(1 * sizeof(int));

    parallel_search_arg->found[0] = 0;
    parallel_search_arg->res_i[0] = -1;

    // Pozivanje radne funkcije za zadatak 1.
    function = working_function_for_task_1;
    argument = parallel_search_arg;

    // Budjenje dretvi
    synchronize(synchronize_all_vars);
    // I glavna dretva izvrsava istu funkciju.
    working_function_for_task_1(0, argument);
    // Sinkroniziranje dretvi
    synchronize(synchronize_all_vars);
    
    // Prekidanje rada dretvi;
    stopped = 1;
    synchronize(synchronize_all_vars);
    // Glavna dretva ceka da kreirane dretve zavrse sa radom.
    for (int i = 0; i < p; i++ ) {
        pthread_join(threads[i], 0);
    }
    
    printf("res_i: %d\n", *parallel_search_arg->res_i);

    for (int i = 0; i < p; i++) {
        free(wd[i]);
    }

    free(wd);
    free(threads);
    free(parallel_search_arg->x);
    free(parallel_search_arg->lim);
    free(parallel_search_arg->c);
    free(parallel_search_arg->q);
    free(parallel_search_arg->found);
    free(parallel_search_arg->res_i);
    free(parallel_search_arg);
    free(synchronize_all_vars);

    return 0;
}