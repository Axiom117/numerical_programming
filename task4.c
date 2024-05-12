#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define FIRST_LINE 11L
#define COLUMN_NUM 5
#define INITIAL 100
#define INTERVAL 5
#define MAXOMEGA 25
#define FILE_NAME4 "task4.csv"
#define FIRST_LINE_OUT4 "threshold,points\n"

/* define the type of data points */
typedef struct {
    double rho;
    double u;
    double v;
    double x;
    double y;
} point_t;

/* file manipulation functions */
FILE* safe_fopen(const char* path, const char* mode)
{
	FILE* fp = fopen(path, mode);
	if (fp == NULL) {
		perror("file open error.");
		exit(EXIT_FAILURE);
	}
	return fp;
}

void exit_if_null(void *ptr, char *msg) {
    if (!ptr) {
        printf("unexpected null pointer: %s\n", msg);
        exit(EXIT_FAILURE);
    }
}

/* this function opens the file can read data */
point_t *read_points(const char *file_name, int *validPoints) {
    size_t currentSize = INITIAL; // set the initial size of the array
    *validPoints = 0;

    point_t *array;
    point_t p; //buffer
    /* initial allocation */
    array = (point_t *)malloc(currentSize * sizeof(point_t));
    FILE *fp = safe_fopen(file_name, "r");
    /* skip the first line */
    fseek(fp, FIRST_LINE, SEEK_SET);

    while (fscanf(fp, "%lf,%lf,%lf,%lf,%lf\n", &p.rho, &p.u, &p.v, &p.x, &p.y) == COLUMN_NUM)
    {
        /* if the space of the array run out */
        if (*validPoints == currentSize) {
            currentSize *= 2;
            /* expand the size of array by two times */
            array = realloc(array, currentSize * sizeof(point_t));
            exit_if_null(array, "reallocation");
        }
        array[*validPoints] = p;
        (*validPoints)++;
    }
    return array;
}

int cmpOmega(const void *v1, const void *v2) {
    if (*(double*)v1 <= *(double*)v2) {
        return 1;
    } else {
        return -1;
    }
}

int main(int argc, char *argv[]) {
    const char *file_name = argv[1];
    point_t **pointGrid;
    point_t *pointArray;
    int distinctX=0, distinctY, validPoints, num = 0;
    double firstY;
    double *omega;

    pointArray = read_points(file_name, &validPoints);

    /* before converting to 2d array, identify how many columns and rows are there in the array */
    firstY = pointArray[0].y;
    for (int i = 0; i < validPoints; i++) {
        /* since we assume the data is sorted in ascending order in y followed by ascending order in x */
        /* we can simply focus on the points that have the same y coordinate */
        if (pointArray[i].y == firstY) {
            distinctX++;
        } else {
            break;
        }
    }
    /* total number of points = m * n */
    distinctY = validPoints / distinctX;

    /* now we have obtained the exact number of columns and rows */
    pointGrid = (point_t **)malloc(distinctX * sizeof(point_t *));
    for (int i=0; i < distinctX; i++) {
        /* allocate space to the secondary structure of pointArray */
        pointGrid[i] = (point_t *)malloc(distinctY * sizeof(point_t));
        for (int j=0; j < distinctY; j++) {
            /* a bit tricky here, because the y is ordered followed by x in pointArray */
            /* then use an equation: distinctX*j + i to jump between the points we want to write in */  
            pointGrid[i][j] = pointArray[distinctX * j + i ];
        }
    }

    free(pointArray);
    pointArray = NULL;

    /* initialize the array w for storing vorticity values */
    omega = (double *)malloc((distinctX - 1) * (distinctY - 1) * sizeof(double));
    /* after finishing the construction of the grid, let's calculate the value of w now */
    for (int i = 0; i < distinctX - 2; i++) {
        for (int j = 0; j < distinctY - 2; j++) {
            omega[num] = (pointGrid[i + 1][j].v - pointGrid[i][j].v) / (pointGrid[i + 1][j].x - pointGrid[i][j].x) - (pointGrid[i][j + 1].u - pointGrid[i][j].u) / (pointGrid[i][j + 1].y - pointGrid[i][j].y);
            num++;
        }
    }
    /* sort the omega array in descending order to make things easier */
    qsort(omega, num, sizeof(double), cmpOmega);
    double threshold;
    int flag=0, n; // in video editing, we usually use "in" and "out" identifiers to crop video
    int *pointsNumber; // record the number of values that are within the threshold range

    n = MAXOMEGA / INTERVAL; // how many samples do we wanna get
    pointsNumber = (int *)malloc(n * sizeof(int));
    threshold = (double)MAXOMEGA;

    pointsNumber[n-1] = 0;
    for (int i = 0; (i < num) && threshold>0; i++) {
        if ((threshold>omega[i]) && (omega[i]>=(threshold-INTERVAL)) && !flag) {
            pointsNumber[n-1]++;
            /* flag = 1 means we've entered the range within the threshold domain */
        } else if (omega[i]<(threshold-INTERVAL)) {
            threshold -= INTERVAL;
            /* we are now out of threshold domain, lower the threshold and reset the flag */
            n -= 1;
            /* initial next number */
            pointsNumber[n - 1] = 0;
        }
    }

    FILE *fp2 = safe_fopen(FILE_NAME4, "w");
    fprintf(fp2, FIRST_LINE_OUT4);
    for (int i = 0; i < MAXOMEGA / INTERVAL; i++) {
        fprintf(fp2, "%d,%d\n", (int)(threshold+INTERVAL), pointsNumber[i]);
        threshold += INTERVAL;
    }
    return 0;
}