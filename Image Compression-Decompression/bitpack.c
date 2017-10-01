/* bitpack.c
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * Bitpack.c packs signed and unsigned uint64_t in a 64 bit word
 * note: homework server uses arithmetic right shift
 */
#include <bitpack.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <except.h>
#include <math.h>

const int MAX_UBITS         = sizeof(uint64_t) * 8;
const int MAX_SBITS         = sizeof(int64_t) * 8;
const int SIZE32            = sizeof(int) * 8;

static inline uint64_t logical_shift_right(uint64_t word, int n);
static inline uint64_t shift_left(uint64_t word, int n);
static inline uint64_t arith_shift_right(uint64_t word, int n);

Except_T Bitpack_Overflow = { "Overflow packing bits" };

/* checks if the given unsigned int would fit in width bits */
bool Bitpack_fitsu(uint64_t n, unsigned width)
{
        /* special case: if width is 0, value has to be 0 */
        if (width == 0) {
                if (n == 0)
                        return 1;
                else
                        return 0;
        }
        /* idea is that values 0 and 1 can fit in anything */
        if (n == (uint64_t)0 || n == (uint64_t)1) {
                return 1;
        }
        unsigned bits;
        /* value can  fit in log2(value) + 1 bits */
        bits = (int) log2(n) + 1;

        if (bits <= width)
                return 1;
        else
                return 0;
}

bool Bitpack_fitss(int64_t n, unsigned width)
{
        /* special case; if width is 0 or 1 (accounting for signed bit), the
         * value has to be 1 */
        if (width == 0 || width == 1) {
                if (n == 0)
                        return 1;
                else
                        return 0;
        }
        /* if value is 0, it can fit in anything */
        if (n == 0) {
                return 1;
        }
        /* value can fit in log(2) + 1 + 1 bits (add one extra to account for
         * signed */
        unsigned bits;
        bits = (int) log2(abs(n)) + 2;

        if (bits <= width)
                return 1;
        else
                return 0;
}

/* get unsigned of width and lsb from word */
uint64_t Bitpack_getu(uint64_t word, unsigned width, unsigned lsb)
{
        /* make a mask that has 1's where number is wanted, then & with word
         * then shift right to obtain value */
        if (width == 0)
                return 0;
        uint64_t mask = ~0;
        mask = shift_left(logical_shift_right(mask, (MAX_UBITS - width)), lsb);
        /* use logical shift here because working with unsigned */
        return logical_shift_right(mask & word, lsb);
}

int64_t Bitpack_gets(uint64_t word, unsigned width, unsigned lsb)
{
        /* make a mask that has 1's where number is wanted, then & with word
         * then shift right to obtain value */
        if (width == 0)
                return 0;
        uint64_t mask = ~0;
        mask = shift_left(logical_shift_right(mask, (MAX_UBITS - width)), lsb);
        /* use arithmetic shift here because working with signed - we want to
         * carry over msb */
        return arith_shift_right(shift_left(mask & word , MAX_UBITS-width-lsb),
                                 MAX_UBITS - width);
}

uint64_t Bitpack_newu(uint64_t word, unsigned width, unsigned lsb,
                      uint64_t value)
{
        /* make a mask over position of desired location by &ing two masks,
         * then & with word and shift value over and | them together to insert
         * value */
        if (!Bitpack_fitsu(value, width)) {
                RAISE(Bitpack_Overflow);
        }
        uint64_t mask;
        mask = ~0;
        mask = shift_left(mask, lsb + width);
        uint64_t mask2 = ~0;
        mask2 = mask2 << lsb;
        mask = mask | ~mask2;
        word = word & mask;
        return (word | (shift_left(value, lsb)));
}
uint64_t Bitpack_news(uint64_t word, unsigned width, unsigned lsb,
                      int64_t value)
{
        /* make a mask over position of desired location by &ing two masks,
         * then & with word and shift value over and | them together to insert
         * value */
        if (!Bitpack_fitss(value, width)) {
                RAISE(Bitpack_Overflow);
        }
        uint64_t ones = ~0;
        ones = shift_left(ones, width);
        uint64_t mask;
        mask = ~0;
        mask = shift_left(mask, lsb + width);
        uint64_t mask2 = ~0;
        mask2 = mask2 << lsb;
        mask = mask | ~mask2;
        word = word & mask;
        /* subtract by ones to remove trailing ones from the value before
         * inserting */
        return (word | (shift_left(value - ones, lsb)));
}

static inline uint64_t logical_shift_right(uint64_t word, int n)
{
        if (n >= MAX_UBITS) {
                //fprintf(stderr, "n %d is greater than %d\n", n, MAX_UBITS);
                return 0;
        }
        else {
                return word >> n;
        }

}

/* logically or arithmetically shifts the given word left
 * returns 0 if shifted 64 bits
 */
static inline uint64_t shift_left(uint64_t word, int n)
{
        if (n >= MAX_UBITS) {
                return 0;
        }
        else {
                return word << n;
        }
}

static inline uint64_t arith_shift_right(uint64_t word, int n)
{
        if (n >= MAX_UBITS) {
                return 0;
        }
        else {
                return ((int64_t) word >> n);
        }
}

