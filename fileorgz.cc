// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "fileorgz.h"

template<typename IO_TYPE>
void FileOrGZ<IO_TYPE>::alloc_buf() {
  buf = (char *) malloc(INIT_SIZE);
  if (buf == NULL) {
    fprintf(stderr, "ERROR: out of memory\n");
    exit(1);
  }
  buf_size = INIT_SIZE;
  buf_len = 0;
}

// open <filename> using standard FILE *
template<>
bool FileOrGZ<FILE *>::open(const char *filename, const char *mode) {
  // First allocate a buffer for I/O:
  alloc_buf();

  fp = fopen(filename, mode);
  if (!fp)
    return false;
  else
    return true;
}

// open <filename> as a gzipped file
template<>
bool FileOrGZ<gzFile>::open(const char *filename, const char *mode) {
  // First allocate a buffer for I/O:
  alloc_buf();

  fp = gzopen(filename, mode);
  if (!fp)
    return false;
  else
    return true;
}

template<>
int FileOrGZ<FILE *>::getline() {
  return ::getline(&buf, &buf_size, fp);
}

template<>
int FileOrGZ<gzFile>::getline() {
  int n_read = 0;
  int c;

  while ((c = gzgetc(fp)) != EOF) {
    // About to have read one more, so n_read + 1 needs to be less than *n.
    // Note that we use >= not > since we need one more space for '\0'
    if (n_read + 1 >= (int) buf_size) {
      const size_t GROW = 1024;
      char *tmp_buf = (char *) realloc(buf, buf_size + GROW);
      if (tmp_buf == NULL) {
	fprintf(stderr, "ERROR: out of memory!\n");
	exit(1);
      }
      buf_size += GROW;
      buf = tmp_buf;
    }
    buf[n_read] = (char) c;
    n_read++;
    if (c == '\n')
      break;
  }

  if (c == EOF && n_read == 0)
    return -1;

  buf[n_read] = '\0';

  return n_read;
}

template<>
int FileOrGZ<FILE *>::printf(const char *format, ...) {
  va_list args;
  int ret;
  va_start(args, format);

  ret = vfprintf(fp, format, args);

  va_end(args);
  return ret;
}

template<>
int FileOrGZ<gzFile>::printf(const char *format, ...) {
  va_list args;
  int ret;
  va_start(args, format);

  // NOTE: one can get automatic parallelization (in a second thread) for
  // gzipped output by opening a pipe to gzip (or bgzip). For example:
  //FILE *pipe = popen("gzip > output.vcf.gz", "w");
  // Can then fprintf(pipe, ...) as if it were a normal file.

  // gzvprintf() is slower than the code below that buffers the output.
  // Saw 13.4% speedup for processing a truncated VCF with ~50k lines and
  // 8955 samples.
//  ret = gzvprintf(fp, format, args);
  ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
  if (ret < 0) {
    printf("ERROR: could not print\n");
    perror("printf");
    exit(10);
  }

  if (buf_len + ret > buf_size - 1) {
    // didn't fit the text in buf
    // first print what was in buf before the vsnprintf() call:
    gzwrite(fp, buf, buf_len);
    buf_len = 0;
    // now ensure that redoing vsnprintf() will fit in buf:
    if ((size_t) ret > buf_size - 1) {
      do { // find the buffer size that fits the last vsnprintf() call
	buf_size += INIT_SIZE;
      } while ((size_t) ret > buf_size - 1);
      free(buf);
      buf = (char *) malloc(buf_size);
      if (buf == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
    }
    // redo:
    ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
  }

  buf_len += ret;
  if (buf_len >= buf_size - 1024) { // within a tolerance of MAX_BUF?
    // flush:
    gzwrite(fp, buf, buf_len);
    buf_len = 0;
  }

  va_end(args);
  return ret;
}

template<>
int FileOrGZ<FILE *>::close() {
  assert(buf_len == 0);
  // should free buf, but I know the program is about to end, so won't
  return fclose(fp);
}

template<>
int FileOrGZ<gzFile>::close() {
  if (buf_len > 0)
    gzwrite(fp, buf, buf_len);
  // should free buf, but I know the program is about to end, so won't
  return gzclose(fp);
}
