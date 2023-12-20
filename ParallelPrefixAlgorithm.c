#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
int counter = 0;
int limit;

void synchronize() {
    pthread_mutex_lock(&mutex);
    counter++;
    if (counter == limit) {
        counter = 0;
        pthread_cond_broadcast(&cond);
    }
    else {
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
}

typedef struct thread_structure {
    int MYPROC;
    int num_threads;
    int n;
    int* x;
    int* y;
    int* local_y;
} ThreadData;


void* parallel_prefix_scan(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    int MYPROC = data->MYPROC;
    int num_threads = data->num_threads;
    int n = data->n;
    int* x = data->x;
    int* y = data->y;
    int* local_y = data->local_y;

    // Grupiranje n elemenata u blokove
    int small = n / num_threads;
    int large = small + 1;
    int num_large = n % num_threads;
    int block = (MYPROC < num_large ? large : small);

    int first, last;
    if (MYPROC < num_large) {
        first = MYPROC * large;
    } 
    else {
        first = num_large * large + (MYPROC - num_large) * small;
    }
    last = first + block;

    // Penjanje
    int k = ceil(log2(n));
    int m = num_threads;
    int step = 1;

    local_y[MYPROC] = 0;
    for (int i = first; i < last; ++i) {
        local_y[MYPROC] += x[i];
    }
    synchronize();    

    for (int j = 2; j <= k; ++j) {
        m = m / 2;

        for (int i = 1; i <= m; ++i) {
            if (MYPROC == i * step - 1) {
                local_y[2 * MYPROC + 1] = local_y[2 * MYPROC - step + 1] + local_y[2 * MYPROC + 1];
            }
        }
        synchronize();

        if (j != k && j != k - 1) step = step * 2;
    }

    // Spuštanje
    m = 1;
    step = step / 2;
    for (int j = 1; j <= k - 2; ++j) {
        m = 2 * m;
        step = step / 2;
        for (int i = 2; i <= m; ++i) {
            if (MYPROC == i * step - 1) {
                local_y[2 * MYPROC - step + 1] = local_y[2 * MYPROC - 2 * step + 1] + local_y[2 * MYPROC - step + 1];
            }
        }
        synchronize();
    }

    k = ceil(log2(block)); 
    int block2 = pow(2, k); 

    int* local_y2 = (int*)malloc(block2 * sizeof(int)); 

    for (int i = 0; i < block; ++i) {
        local_y2[i] = x[first + i];
    }

    for (int i = block; i < block2; ++i) {
        local_y2[i] = 0;
    }

    int sum = 0;
    for (int i = 0; i < block2; i++) {
        sum += local_y2[i];
        local_y2[i] = sum;
    }
    synchronize();

    // Spremanje rezultata u y, osim za nultu dretvu
    for (int i = first; i < last; ++i) {
        y[i] = local_y2[i - first] + (i > 0 ? local_y[MYPROC - 1] : 0);
    }
    synchronize();

    free(local_y2); 

    return NULL;
}

int main() {
    int num_threads = 16;
	int n = 100;
	pthread_t *threads;
	ThreadData *arg;
	int *x, *y, *local_y;

	limit = num_threads;
	threads = (pthread_t*) malloc((num_threads)*sizeof(pthread_t));
	x = (int*) malloc(n*sizeof(int));
	y = (int*) malloc(n*sizeof(int));
	local_y = (int*) malloc(num_threads*sizeof(int));
	arg = (ThreadData*) malloc(num_threads*sizeof(ThreadData));

    // Inicijalizacija ulaznog vektora x (samo jedinice)
    for (int i = 0; i < n; ++i) {
        x[i] = 1;
    }

    // Stvaranje dretvi
    for (int i = 0; i < num_threads - 1; ++i) {
		(arg+i+1)->num_threads = num_threads;
		(arg+i+1)->MYPROC = i + 1;
		(arg+i+1)->n = n;
		(arg+i+1)->x = x;
		(arg+i+1)->y = y;
		(arg+i+1)->local_y = local_y;

        pthread_create( &threads[i], NULL, parallel_prefix_scan, (void*) (arg+i+1));
    }

	arg->num_threads = num_threads;
	arg->MYPROC = 0;
	arg->n = n;
	arg->x = x;
	arg->y = y;
	arg->local_y = local_y;
	parallel_prefix_scan( (void*) arg );

    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

    // Ispis rezultata lokalnih parcijalnih suma
    printf("Izlazni vektor local_y:\n");
    for (int i = 0; i < num_threads; ++i) {
        printf("%d ", local_y[i]);
    }
    printf("\n");

    // Ispis rezultata vektora y
    printf("Izlazni vektor y:\n");
    for (int i = 0; i < n; ++i) {
        printf("%d ", y[i]);
    }
    printf("\n");

    free(threads);
	free(x);
	free(y);
	free(local_y);
	free(arg);
	
    return 0;
}