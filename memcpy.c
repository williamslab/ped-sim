#include <string.h>

// some systems do not have newest memcpy@@GLIBC_2.14 - stay with old good one
// see https://stackoverflow.com/questions/8823267/linking-against-older-symbol-version-in-a-so-file
ifdef __linux__
asm (".symver memcpy, memcpy@GLIBC_2.2.5");
#endif
void *__wrap_memcpy(void *dest, const void *src, size_t n)
{
    return memcpy(dest, src, n);
}
