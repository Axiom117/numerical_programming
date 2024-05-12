#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define FIRST_LINE 11L
#define INITIAL 100 // let's say, set the initial size of the array 100
#define COLUMN_NUM 5
#define CENTRELINE 0
#define PERCENTAGE 0.4
#define FILE_NAME3 "task3.csv"


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

/* define the balance factor of BST */
typedef enum {
    EH = 0,
    LH = 1,
    RH = -1
} bh_t;

/* define the BST node, however, I can't find a better name, I simply name it "leaf" */
typedef struct leaf leaf_t;
struct leaf {
    point_t point;
    int bf; // balance factor
    leaf_t *left, *right;
};

/* define the BST */
typedef struct {
    leaf_t *root;
    int (*cmp)(void *, void *);
    void (*del)(void *);
} tree_t;

/* some BST manipulation functions */
tree_t *make_tree(void fdel(void*)) {
    tree_t *tree;
    tree = malloc(sizeof(*tree));
    /* initialize the tree */
    tree->root = NULL;
    /* refer to the function pointer */
    tree->del = fdel;
    return tree;
}

/* function for delete single tree node */
void delTree(void *v) {
    leaf_t *leaf = (leaf_t *)v;
    if (leaf == NULL) {
        return;
    }
    delTree(leaf->left);
    delTree(leaf->right);
    free(leaf);
}

/* the functions that allocate array into BST are super-complex, however, when it comes to sorted array, things are getting much easier */ 
leaf_t *array_to_BST(point_t* array, int start, int end) {
    if (start > end) {
        return NULL;
    }
    int mid = start + (end - start) / 2;
    leaf_t *leaf = (leaf_t *)malloc(sizeof(leaf_t *));
    leaf->point = array[mid];
    /* recursively construct the binary structure */
    leaf->left = array_to_BST(array, start, mid - 1);
    leaf->right = array_to_BST(array, mid + 1, end);
    return leaf;
}

leaf_t *tree_search(tree_t *t, leaf_t* leaf, point_t p) {
    if(leaf) {
        /* go left as p < leaf->point */
        if (t->cmp(&p, &leaf->point) == -1) {
            return tree_search(t, leaf->left, p);
        }
        /* go right as p > leaf->point */
        if (t->cmp(&p, &leaf->point) == 1) {
            return tree_search(t, leaf->right, p);
        }
    }
    return leaf;
}

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



int main(int argc, char *argv[]) {
    const char *file_name = argv[1];
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
    pointTree = make_tree(delTree);

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
    pointTree->root = array_to_BST(pointArray, 0, validPoints - 1);

    FILE *fp2 = safe_fopen(FILE_NAME3, "w");
    /* Linear search on the array */
    linear_search_array(fp2, pointArray, target, validPoints);
    fprintf(fp2, "\n");
    /* Binary search on the array */
    binary_search(fp2, pointArray, target, 0, validPoints-1);
    fprintf(fp2, "\n");
    /* Linear search on the Linked list */
    linear_search_list(fp2, pointList, target);
    /* Search on the balanced BST */

    free(pointArray);
    pointArray = NULL;

    /* no need to manully free the list as the elements within has been popped out every time searching */

    fclose(fp);
    return 0;
}