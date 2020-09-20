// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <zlib.h>

#ifndef FILEORGZ_H
#define FILEORGZ_H

template<typename IO_TYPE>
class FileOrGZ {
  public:
    bool open(const char *filename, const char *mode);
    int getline();
    int printf(const char *format, ...);
    int close();

    static const int INIT_SIZE = 1024 * 50;

    // IO_TYPE is either FILE* or gzFile;
    IO_TYPE fp;

    // for I/O:
    char *buf;
    size_t buf_size;
    size_t buf_len;

  private:
    void alloc_buf();
};

#endif // FILEORGZ_H
