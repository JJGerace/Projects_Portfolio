/* type_manip.h
 * by Taher Mun and Jacob Gerace, 10/16/2013
 *
 * type_manip.c converts the 4 parameters a,b,c,d to and from
 * their relevant integer type, can also be reused to work with
 * other conversions (ex: a double to a 10 bit integer instead of
 * a 9 bit) by changing constants
 * */

int64_t double_to_signed(double d);
int64_t double_bcd_to_signed(double d);
uint64_t double_a_to_unsigned(double d);
double signed_bcd_to_double(int64_t i);
double unsigned_a_to_double(uint64_t i);

