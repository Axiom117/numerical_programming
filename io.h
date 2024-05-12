#ifndef __IO_H__
#define __IO_H__
#include <stdio.h> 

   int file_exists(const char* path);
   void* safe_malloc(size_t num_bytes);
   FILE* safe_fopen(const char* path, const char* mode);
   size_t file_size_in_bytes(FILE* fp);
   
#endif