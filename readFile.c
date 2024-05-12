#include <stdio.h>
#include <stdlib.h>

#define LINE_LIMIT 2
#define SEPARATOR "==============================="

void first_lines(FILE *fp, int n);

int main(int argc, char *argv[]) {
    int fnum;
    FILE *fp;

    for (fnum = 1; fnum < argc; fnum++) {
        fprintf(stderr, "Opening %s: ", argv[fnum]);
        if ((fp = fopen(argv[fnum], "r"))==NULL) {
            fprintf(stderr, "........failed\n");
        }else {
            fprintf(stderr, "\n");
            printf("%s %s\n", SEPARATOR, argv[fnum]);
            first_lines(fp, LINE_LIMIT);
            fclose(fp);
        }
    }
    return 0;
}

void
first_lines(FILE *fp, int n) {
    int c;
    int lines = 0;
    while ((c = getc(fp)) != EOF) {
        if (lines < n) {
            putchar(c);
        }
        lines += (c == '\n');
    }
    printf("[%d lines in total]\n", lines);
    return;
}