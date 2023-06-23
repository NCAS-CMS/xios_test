#ifdef USE_UDUNITS
#include <udunits2.h>

ut_system *read_xml(const char *path)
{
  ut_system *utsystem;

  ut_set_error_message_handler(ut_ignore);

  utsystem = ut_read_xml(path);
  if (utsystem == NULL)
  {
    ut_set_error_message_handler(ut_write_to_stderr);
    utsystem = ut_read_xml(path);
    return NULL;
  }

  return utsystem;
}

int are_convertible(ut_system *utsystem, const char *unit1, const char *unit2)
{
  ut_unit *utunit1, *utunit2;
  int ret;

  utunit1 = ut_parse(utsystem, unit1, UT_ASCII);
  utunit2 = ut_parse(utsystem, unit2, UT_ASCII);

  ret = ut_are_convertible(utunit1, utunit2);

  ut_free(utunit1);
  ut_free(utunit2);

  return ret;
}
#endif
