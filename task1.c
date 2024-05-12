#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <io.h>

#define MAX_BUF_LEN 2048
#define XCOORDI 20.00
#define COLUMN_NUM 5
#define POINT_NUM 4
#define FILE_NAME "task1.csv"

FILE* safe_fopen(const char* path, const char* mode)
{
	FILE* fp = fopen(path, mode);
	if (fp == NULL) {
		perror("file open error.");
		exit(EXIT_FAILURE);
	}
	return fp;
}

/* define the type of data points */
typedef struct {
    double rho;
    double u;
    double v;
    double x;
    double y;
} point_t;

int main(int argc, char *argv[]) {
    const char* file_name = argv[1];
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

	FILE* fp2 = safe_fopen(FILE_NAME, "w"); // this creates a new file where the data will be written into
    
	/* write the data into the new file */
	for (int i = 0; i < POINT_NUM; i++) {
		fprintf(fp2, "%lf,%lf,%lf,%lf,%lf\n", points[i].rho, points[i].u, points[i].v,
             points[i].x, points[i].y);
    }

	fclose(fp);
	fclose(fp2);

    return 0;
}