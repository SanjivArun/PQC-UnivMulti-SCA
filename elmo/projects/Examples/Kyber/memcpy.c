/* Minimal memcpy for ELMO - uses no-stack optimization */
#include <stdint.h>
__attribute__((optimize("O2", "omit-frame-pointer", "no-stack-protector")))
void *memcpy(void *dst, const void *src, uint32_t n) {
    const uint8_t *s = (const uint8_t *)src;
    uint8_t *d = (uint8_t *)dst;
    uint32_t i;
    for (i = 0; i < n; i++)
        d[i] = s[i];
    return (void *)d;
}
