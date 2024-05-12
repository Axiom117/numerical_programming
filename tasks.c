/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 926204
 *   Name        : HAORAN YAO
 *
 ***************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "tasks.h"
#include <stdint.h>
#include <assert.h>
#include <sys/time.h>

#define MAX_BUF_LEN 2048
#define XCOORDI 20.00
#define COLUMN_NUM 5
#define POINT_NUM 4
#define FILE_NAME1 "task1.csv"
#define ZERO 0
#define DECIMAL 1000000
#define LEFT_BORDER -15
#define RIGHT_BORDER 85
#define UPPER_BORDER 20
#define LOWER_BORDER -20
#define INITIAL 100 //the initial size of a single row in rows array
#define FIRST_LINE 11L // skip the first line of csv file
#define FIRST_LINE_OUT "rho, u, v, x, y, S\n" //the first line in output file
#define FILE_NAME2 "task2.csv"
#define INITIAL 100 // let's say, set the initial size of the array 100
#define CENTRELINE 0
#define PERCENTAGE 0.4
#define FILE_NAME3 "task3.csv"
#define INTERVAL 5
#define MAXOMEGA 25
#define FILE_NAME4 "task4.csv"
#define FIRST_LINE_OUT4 "threshold,points\n"
#define MILLI2MICRO 1000.0 // convert millisecond to microsecond

/* define the type of data points */
typedef struct {
    double rho;
    double u;
    double v;
    double x;
    double y;
} point_t;

/* now, let declare and define the linked list and node */
typedef struct node node_t;
struct node {
    point_t point;
    node_t *next;
};
typedef struct {
    node_t *head;
    node_t *foot;
} point_list;

/* define the type of cells */
typedef struct {
    point_list *points;
    int k;
    double avRho;
    double avU;
    double avV;
    double avX;
    double avY;
    double score;
} cell_t;

/* define the BST node, however, I can't find a better name, I simply name it "leaf" */
typedef struct tree tree_t;
struct tree {
    point_t point;
    tree_t *left, *right;
};

/* function for delete single tree node */
void delTree(tree_t *root) {
    if (root == NULL) {
        return;
    }
    delTree(root->left);
    delTree(root->right);
    free(root);
    root = NULL;
}

/* the functions that allocate array into BST are super-complex, however, when it comes to sorted array, things are getting much easier */ 
tree_t *array_to_BST(point_t* array, int start, int end) {
    if (start > end) {
        return NULL;
    }
    tree_t *new;
    new = (tree_t*)malloc(sizeof(*new));
    int mid = start + (end - start) / 2;
    new->point = array[mid];
    /* recursively construct the binary structure */
    new->left = array_to_BST(array, start, mid - 1);
    new->right = array_to_BST(array, mid + 1, end);
    return new;
}

void tree_search(FILE *fp, tree_t* leaf, double target, double minDiff, double *result) {
    if(leaf) {
        double diff, flux = (leaf->point.u) * (leaf->point.rho);
        fprintf(fp, "%lf, ", flux);
        diff = fabs(flux - target);
        /* track the minimum difference and record the value of flux */
        if (diff < minDiff) {
            minDiff = diff;
            *result = flux;
        }
        /* go left as target < leaf->point.flux */
        if (flux>target) {
            return tree_search(fp, leaf->left, target, minDiff, result);
        }
        /* go right as target > leaf->point.flux */
        if (flux<target) {
            return tree_search(fp, leaf->right, target, minDiff, result);
        }
    }
}

/* set of file manipulation functions */
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

/* several linked list manipulating functions */
point_list *make_empty_list (void) {
    point_list *list;
    list = (point_list *)malloc(sizeof(*list));
    list -> head = list -> foot = NULL;
    return list;
}

int is_empty_list (point_list *list) {
    return list -> head == NULL;
}

void free_list (point_list *list) {
    node_t *curr, *prev;
    curr = list->head;
    while (curr) {
      prev = curr;
      curr = curr->next;
      free(prev);
    }
    free(list);
}

point_list *insert_at_head(point_list *list, point_t value) {
    node_t *new;
    new = (node_t*)malloc(sizeof(*new));
    new->point = value;
    new->next = list->head;
    list->head = new;
    if (list->foot==NULL) {
        /* if this is the first insertion into the list */
        list->foot = new;
    }
    return list;
}

point_list *insert_at_foot(point_list *list, point_t value) {
    node_t *new;
    new = (node_t*)malloc(sizeof(*new));
    new->point = value;
    new->next = list->head;
    if (list->foot==NULL) {
        /* if this is the first insertion into the list */
        list->head = list->foot = new;
    } else {
        list->foot->next = new;
        list->foot = new;
    }
    return list;
}

point_t get_head(point_list *list) {
    return list->head->point;
}

point_list *get_tail(point_list *list) {
    node_t *oldhead;
    oldhead = list->head;
    list->head = list->head->next;
    if (list->head==NULL) {
        /* if there only one node in the list */
        list->foot = NULL;
    }
    free(oldhead);
    return list;
}

/* sets of functions of quick sort */

int cmpScore(cell_t *c1, cell_t *c2) {
    if ((*c1).score <= (*c2).score) {
        return 1;
    } else {
        return -1;
    }
}

void swap_cell(cell_t *c1, cell_t *c2) {
    cell_t temp;
    temp = *c1;
    *c1 = *c2;
    *c2 = temp;
}

void partition(cell_t A[], int n, cell_t *pivot, int *first_eq, int *first_gt) {
    int next = 0, fe = 0, fg = n, outcome;
    while (next<fg) {
        if ((outcome = cmpScore(A+next, pivot)) < 0) {
            swap_cell(A+fe, A+next);
            fe += 1;
            next += 1;
        } else if (outcome > 0) {
            fg -= 1;
            swap_cell(A + next, A + fg);
        } else {
            next += 1;
        }
    }
    *first_eq = fe;
    *first_gt = fg;
    return;
}

cell_t choose_pivot(cell_t A[], int n) {
    /* consider the cost of evaluating the random number generator is relatively high, we simple pick the middle data */
    return A[n / 2];
}
void quick_sort(cell_t A[], int n) {
    cell_t pivot;
    int first_eq, first_gt;
    if (n<=1) {
        return;
    }
    pivot = choose_pivot(A, n);
    partition(A, n, &pivot, &first_eq, &first_gt);
    quick_sort(A, first_eq);
    quick_sort(A + first_gt, n - first_gt);
}

/* this function specifies the comparison criteria, it will thereby be passed as a function pointer in both qsort() function */
int cmpFlux(const void *v1, const void *v2) {
    if (((*(point_t*)v1).rho*(*(point_t*)v1).u) >= ((*(point_t*)v2).rho*(*(point_t*)v2).u)) {
        return 1;
    } else {
        return -1;
    }
}

/* the main idea of finding the closest value is using the advantage of sorting */
/* since the array (linked list) is sorted in ascending order already */
/* start searching from the smallest value, the absolute difference between the current value and target value must be descending until reach the closest value, then the difference would ascend */
/* we only need to catch this turning point where the minimum absolute difference occurs */
void linear_search_array(FILE *fp, point_t* p, double target, int length) {
    double curDiff, nextDiff, flux;

    for (int i = 0; i < length-1; i++) {
        flux = p[i].u * p[i].rho;
        curDiff = fabs(target - p[i].u * p[i].rho);
        nextDiff = fabs(target - p[i+1].u * p[i+1].rho);
        /* catch the turning point */
        if (curDiff < nextDiff) {
            fprintf(fp, "%lf", flux);
            break;
        }
        fprintf(fp, "%lf, ", flux);
    }
}

/* for linked list, things are getting a bit complicated */
/* since we can't foresee the data of the next point as we did in search_array */
/* we will focus on previous point */
void linear_search_list(FILE *fp, point_list *list, double target) {
    double curDiff, preDiff, flux, preFlux;
    int flag = 0;
    point_t p;
    while (!is_empty_list(list)) {
        p = get_head(list);
        flux = p.u * p.rho;
        curDiff = fabs(target - flux);
        if (!flag) {
            preDiff = curDiff; // initialize the previous difference
            flag = 1;
            preFlux = flux;
            list = get_tail(list);
            continue;
        } else {
            fprintf(fp, "%lf", preFlux);
        }
        if (curDiff > preDiff) {
            break;
        }
        fprintf(fp, ", ");
        preFlux = flux;
        preDiff = curDiff;
        list = get_tail(list);
    }
}

/* main idea for Binary Search is recognizing the "slope" of differences curves between the current flux and the target */
/* aforementioned feature of this sorted array is the concave difference curve */
/* the differences first go down, then rise */

int binary_search(FILE *fp, point_t * p, double target, double left, double right)
{
    int midIndex = (left + right + 1) / 2;
    /* pick three points here to ensure correct comparison */
    double curDiff = fabs(target - p[midIndex].u * p[midIndex].rho);
    double preDiff = fabs(target - p[midIndex-1].u * p[midIndex-1].rho);
    double nextDiff = fabs(target - p[midIndex+1].u * p[midIndex+1].rho);

    double flux = p[midIndex].rho * p[midIndex].u;
    /* if we are at the right place */
    if ((curDiff-preDiff)<0.0 && (nextDiff-curDiff)>0.0) {
        fprintf(fp, "%lf", flux);
        return 1;
    /* therefore, if the gradient is negative, it means we are on the left of what we are searching for */
    } else if ((curDiff-preDiff)<0.0 && (nextDiff-curDiff)<0.0) {
        fprintf(fp, "%lf, ", flux);
        return binary_search(fp, p, target, midIndex, right);
    /* if the gradient is positive, we are on the right */
    } else {
        fprintf(fp, "%lf, ", flux);
        return binary_search(fp, p, target, left, midIndex);
    }
}

/* this function opens the file and read data */
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
    fclose(fp);
    return array;
}

int cmpOmega(const void *v1, const void *v2) {
    if (*(double*)v1 <= *(double*)v2) {
        return 1;
    } else {
        return -1;
    }
}

void maxfluxdiff(const char* flow_file)
{
    const char* file_name = flow_file;
	point_t p; //act as a buffer
	double uflux, vflux, uMax=0, uMin, vMax=0, vMin;
	int flag = 1; //help to initialize the minimum flux
    point_t points[4]; //record the data of points that satisifies the maximum difference
	
    FILE* fp = safe_fopen(file_name, "r");

	/* skip the first row */
    fseek(fp, 11L, SEEK_SET);

    while (fscanf(fp, "%lf,%lf,%lf,%lf,%lf\n", &p.rho, &p.u, &p.v, &p.x, &p.y) == COLUMN_NUM)
    {   
		/* the initialization of min will only be executed once */
		if (flag) {
			uMin = p.rho*p.u;
			vMin = p.rho*p.v;
			flag = 0;
		}
		/* check if x is greater than 20 */
        if (p.x > XCOORDI) {
			uflux = p.rho*p.u;
			vflux = p.rho*p.v;
			if (uflux > uMax) {
				points[0] = p;
				uMax = uflux;
            } else if (uflux < uMin) {
                points[1] = p;
                uMin = uflux;
            } 
			if (vflux > vMax) {
				points[2] = p;
				vMax = vflux;
			} else if (vflux < vMin) {
				points[3] = p;
				vMin = vflux;
			}
        }
    }

	FILE* fp2 = safe_fopen(FILE_NAME1, "w"); // this creates a new file where the data will be written into
    
	/* write the data into the new file */
	for (int i = 0; i < POINT_NUM; i++) {
		fprintf(fp2, "%lf,%lf,%lf,%lf,%lf\n", points[i].rho, points[i].u, points[i].v,
             points[i].x, points[i].y);
    }

	fclose(fp);
	fclose(fp2);
}

void coarsegrid(const char* flow_file, int resolution)
{
    const char *file_name = flow_file;
    point_t o, p, q;   //buffers
    point_list **rows; // create row-arrays for storing points by their y-coordinates from up to down
    int r, c, cActual; // location judging indeices, r plays as the row number, while c stans for the cell number, cActual is a tricky one, I will describe it later
    cell_t *cells; // we want all cells to form an one-dimensional array.

    const int xDivision = (RIGHT_BORDER-LEFT_BORDER)/resolution;
    const int yDivision = (UPPER_BORDER-LOWER_BORDER)/resolution;

    /* initialize rows */
    rows = (point_list **)malloc(resolution*sizeof(point_list *));
    for (int i=0; i<resolution; i++) {
        rows[i] = make_empty_list();
    }
    exit_if_null(rows, "initial allocation!");

    /* initialize cells */
    cells = (cell_t *)malloc(resolution*resolution*sizeof(cell_t));
    for (int j=0; j<resolution*resolution; j++) {
        cells[j].points = make_empty_list();
    }

    /* open the csv file */
    FILE* fp = safe_fopen(file_name, "r");
	/* skip the first row */
    fseek(fp, FIRST_LINE, SEEK_SET);

    /* first divide the points into rows by their y coordinates */
    while (fscanf(fp, "%lf,%lf,%lf,%lf,%lf\n", &p.rho, &p.u, &p.v, &p.x, &p.y) == COLUMN_NUM)
    {   
        /* exclude those poiints that lie on the grid line */
        if (((int)(DECIMAL*(UPPER_BORDER-(p.y)))%(DECIMAL*yDivision) == ZERO) 
            || (((int)(DECIMAL*(p.x-LEFT_BORDER))%(DECIMAL*xDivision) == ZERO))) {
            continue; // why using "DECIMAL" here? To make the decimal number completely integer to avoid auto-round-up or down.
        }
        r = ((double)UPPER_BORDER-p.y)/yDivision; // very essential, determine which row the point belongs, also matches the row number in row array
        rows[r] = insert_at_head(rows[r], p);
    }
    /* then focus on each row */
    for (int i=0; i<resolution; i++) {
        /* traverse through the rows and allocate points into cells by their x coordinates */
        while (!is_empty_list(rows[i])) {
          q = get_head(rows[i]);
          /* same, filter the points that are on vertical grid lines */

          c = (q.x - (double)LEFT_BORDER) /
              xDivision; // c matches the column order, but not the actual cell
                         // number
          /* since the cell array is one dimensional, there must be a data
           * convertion, the cActual determines which cell the point ought to be
           * assigned */
          cActual = c + i * resolution; // quiet like carry-in mechanism
          cells[cActual].points = insert_at_head(cells[cActual].points, q);
          rows[i] = get_tail(rows[i]);
        }
        /* free the list to reserve space */
        free_list(rows[i]);
        rows[i] = NULL;
    }
    free(rows);
    rows = NULL;

    /* now do the calculation within each cell */
    int i;
    double totRho, totX, totY, totU, totV;
    for (i=0; i < resolution*resolution; i++) {
        /* intialize cell data */
        cells[i].k = 0;
        totRho = 0;
        totU = 0;
        totV = 0;
        totX = 0;
        totY = 0;
        while (!is_empty_list(cells[i].points)) {
            o = get_head(cells[i].points);
            totRho += o.rho;
            totX += o.x;
            totY += o.y;
            totU += o.u;
            totV += o.v;
            cells[i].k++;
            cells[i].points = get_tail(cells[i].points);
        }
        if (cells[i].k != 0) {
            cells[i].avRho = (totRho / (double)cells[i].k);
            cells[i].avX = (totX/(double)cells[i].k);
            cells[i].avY = (totY/(double)cells[i].k);
            cells[i].avU = (totU/(double)cells[i].k);
            cells[i].avV = (totV/(double)cells[i].k);
            cells[i].score = 100 * (sqrt(cells[i].avU * cells[i].avU + cells[i].avV * cells[i].avV)) / (sqrt(cells[i].avX * cells[i].avX + cells[i].avY * cells[i].avY));
        }
        /* free cells' points data to reserve space, only averages and score are needed */
        free_list(cells[i].points);
        cells[i].points = NULL;
    }

    quick_sort(cells, i);

    FILE* fp2 = safe_fopen(FILE_NAME2, "w");
    fprintf(fp2, FIRST_LINE_OUT);
    for (int j = 0; j < i; j++) {
        fprintf(fp2, "%lf,%lf,%lf,%lf,%lf\n",  cells[j].avU, cells[j].avV, cells[j].avX, cells[j].avY, cells[j].score);
    }

    free(cells);
    cells = NULL;

    fclose(fp);
    fclose(fp2);
}

void searching(const char* flow_file)
{
    struct timeval start;
	struct timeval stop;

    const char *file_name = flow_file;
    size_t currentSize = INITIAL; // set the initial size of the array
    int validPoints = 0;
    point_t *pointArray; // array for storing points
    point_list *pointList; // linked list for storing points
    tree_t *pointTree; // balanced search tree

    point_t p; // buffer
    double target, maxFlux=0.0;

    /* initialize the array */
    pointArray = (point_t *)malloc(currentSize * sizeof(point_t));
    exit_if_null(pointArray, "initial allocation!");
    /* initialize the linked list */
    pointList = make_empty_list();
    /* initialize the BST */
    pointTree = (tree_t *)malloc(sizeof(tree_t*));

    /* now open the file */
    FILE *fp = safe_fopen(file_name, "r");
    /* skip the first row */
    fseek(fp, FIRST_LINE, SEEK_SET);

     while (fscanf(fp, "%lf,%lf,%lf,%lf,%lf\n", &p.rho, &p.u, &p.v, &p.x, &p.y) == COLUMN_NUM)
    {
        /* if the space of the array run out */
        if (validPoints == currentSize) {
            currentSize *= 2;
            /* expand the size of array by two times */
            pointArray = realloc(pointArray, currentSize * sizeof(point_t));
            exit_if_null(pointArray, "reallocation");
        }
        /* check if the point is on the centre line */
        if (p.y == (double)CENTRELINE) {
            pointArray[validPoints] = p;
            validPoints++;
            if (p.rho*p.u > maxFlux) {
                maxFlux = p.rho * p.u;
            }
        }
    }

    /* utilize the polymorphism of qsort, pass the cmpFlux as the funtion pointer to give hugh flexbility of sorting */
    qsort(pointArray, validPoints, sizeof(*pointArray), cmpFlux);

    /* compute the value of 40% of the maximum flux */
    target = maxFlux * PERCENTAGE;

    /* data writing */
    for (int i = 0; i < validPoints; i++) {
        /* insert the sorted points to linked list */
        pointList = insert_at_foot(pointList, pointArray[i]);
    }
    /* write sorted array to BST */
    //pointTree = array_to_BST(pointArray, 0, validPoints - 1);
 
    FILE *fp2 = safe_fopen(FILE_NAME3, "w");

    /* Linear search on the array */
    gettimeofday(&start, NULL);
    linear_search_array(fp2, pointArray, target, validPoints);
    fprintf(fp2, "\n");
    gettimeofday(&stop, NULL);
    double elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
    printf("TASK 3 Array Linear Search:  %.2f microseconds\n", elapsed_ms*MILLI2MICRO);

    /* Binary search on the array */
    gettimeofday(&start, NULL);
    binary_search(fp2, pointArray, target, 0, validPoints-1);
    fprintf(fp2, "\n");
    gettimeofday(&stop, NULL);
    elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
    printf("TASK 3 Array Binary Search:  %.2f microseconds\n", elapsed_ms*MILLI2MICRO);

    /* Linear search on the Linked list */
    gettimeofday(&start, NULL);
    linear_search_list(fp2, pointList, target);
    fprintf(fp2, "\n");
    gettimeofday(&stop, NULL);
    elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
    printf("TASK 3 List Linear Search:  %.2f microseconds\n", elapsed_ms*MILLI2MICRO);
    
    free(pointArray);
    pointArray = NULL;

    /* no need to manully free the list as the elements within has been popped out every time searching */
    free(pointList);

    //delTree(pointTree);

    fclose(fp);
    fclose(fp2);
}

void vortcalc(const char* flow_file)
{
    const char *file_name = flow_file;
    point_t **pointGrid;
    point_t *pointArray;
    int distinctX=0, distinctY, validPoints=0, num = 0;
    double firstY;
    double *omega;

    size_t currentSize = INITIAL; // set the initial size of the array

    point_t p; //buffer
    /* initial allocation */
    pointArray = (point_t *)malloc(currentSize * sizeof(point_t));
    FILE *fp = safe_fopen(file_name, "r");
    /* skip the first line */
    fseek(fp, FIRST_LINE, SEEK_SET);

    while (fscanf(fp, "%lf,%lf,%lf,%lf,%lf\n", &p.rho, &p.u, &p.v, &p.x, &p.y) == COLUMN_NUM)
    {
        /* if the space of the array run out */
        if (validPoints == currentSize) {
            currentSize *= 2;
            /* expand the size of array by two times */
            pointArray = realloc(pointArray, currentSize * sizeof(point_t));
            exit_if_null(pointArray, "reallocation");
        }
        pointArray[validPoints] = p;
        validPoints++;
    }

    //pointArray = read_points(file_name, &validPoints);

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

    /* free the 2d array */
    for (int i = 0; i < distinctX; i++)
    {
        point_t* currentIntPtr = pointGrid[i];
        free(currentIntPtr);
    }
    
    free(pointGrid);
    pointGrid = NULL;

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

    free(omega);

    FILE *fp2 = safe_fopen(FILE_NAME4, "w");
    fprintf(fp2, FIRST_LINE_OUT4);
    for (int i = 0; i < MAXOMEGA / INTERVAL; i++) {
        fprintf(fp2, "%d,%d\n", (int)(threshold+INTERVAL), pointsNumber[i]);
        threshold += INTERVAL;
    }

    free(pointsNumber);

    fclose(fp);
    fclose(fp2);
}
