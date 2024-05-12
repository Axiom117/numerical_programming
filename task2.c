#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define ZERO 0
#define DECIMAL 1000000
#define COLUMN_NUM 5
#define RESOLUTION 10
#define LEFT_BORDER -15
#define RIGHT_BORDER 85
#define UPPER_BORDER 20
#define LOWER_BORDER -20
#define INITIAL 100 //the initial size of a single row in rows array
#define FIRST_LINE 11L // skip the first line of csv file
#define FIRST_LINE_OUT "rho, u, v, x, y, S\n" //the first line in output file
#define FILE_NAME "task2.csv"

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

int main(int argc, char *argv[]) {
    const char* file_name = argv[1];
    point_t o, p, q; //buffers
    point_list **rows; // create row-arrays for storing points by their y-coordinates from up to down
    int r, c, cActual; // location judging indeices, r plays as the row number, while c stans for the cell number, cActual is a tricky one, I will describe it later
    cell_t *cells; // we want all cells to form an one-dimensional array.

    const int xDivision = (RIGHT_BORDER-LEFT_BORDER)/RESOLUTION;
    const int yDivision = (UPPER_BORDER-LOWER_BORDER)/RESOLUTION;

    /* initialize rows */
    rows = (point_list **)malloc(RESOLUTION*sizeof(point_list *));
    for (int i=0; i<RESOLUTION; i++) {
        rows[i] = make_empty_list();
    }
    exit_if_null(rows, "initial allocation!");

    /* initialize cells */
    cells = (cell_t *)malloc(RESOLUTION*RESOLUTION*sizeof(cell_t));
    for (int j=0; j<RESOLUTION*RESOLUTION; j++) {
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
    for (int i=0; i<RESOLUTION; i++) {
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
          cActual = c + i * RESOLUTION; // quiet like carry-in mechanism
          cells[cActual].points = insert_at_head(cells[cActual].points, q);
          rows[i] = get_tail(rows[i]);
        }
        /* free the list to reserve space */
        free_list(rows[i]);
        rows[i] = NULL;
    }

    /* now do the calculation within each cell */
    int i;
    double totRho, totX, totY, totU, totV;
    for (i=0; i < RESOLUTION*RESOLUTION; i++) {
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

    FILE* fp2 = safe_fopen(FILE_NAME, "w");
    fprintf(fp2, FIRST_LINE_OUT);
    for (int j = 0; j < i; j++) {
        fprintf(fp2, "%lf,%lf,%lf,%lf,%lf\n",  cells[j].avU, cells[j].avV, cells[j].avX, cells[j].avY, cells[j].score);
    }

    free(cells);
    cells = NULL;

    fclose(fp);
    fclose(fp2);

    return 0; 
}