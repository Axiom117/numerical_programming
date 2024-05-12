#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>

int file_exists(const char* path)
{
	FILE* fp = fopen(path, "r");
	if (fp == NULL) {
		return 0;
	}
	fclose(fp);
	return 1;
}

void* safe_malloc(size_t num_bytes)
{
	void* ptr = malloc(num_bytes);
	if (ptr == NULL) {
		printf("ERROR: malloc(%zu)\n", num_bytes);
		exit(EXIT_FAILURE);
	}
	return ptr;
}

FILE* safe_fopen(const char* path, const char* mode)
{
	FILE* fp = fopen(path, mode);
	if (fp == NULL) {
		perror("file open error.");
		exit(EXIT_FAILURE);
	}
	return fp;
}

size_t file_size_in_bytes(FILE* fp)
{
	assert(fp != NULL);
	long cur_pos = ftell(fp);
	if (cur_pos == -1) {
		perror("ftell error.");
		exit(EXIT_FAILURE);
	}
	/* seek to the end */
	int ret = fseek(fp, 0L, SEEK_END);
	if (ret == -1) {
		perror("file seek error.");
		exit(EXIT_FAILURE);
	}
	long fs = ftell(fp);

	/* seek back to the previous position */
	ret = fseek(fp, cur_pos, SEEK_SET);
	if (ret == -1) {
		perror("file seek error.");
		exit(EXIT_FAILURE);
	}
	return (size_t)fs;
}