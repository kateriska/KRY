#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <assert.h>
#include <bits/stdc++.h>
#include <list>
#include <vector>
#include <gmpxx.h>
#include <getopt.h>

using namespace std;

bool miller_rabin_test(mpz_t n)
{
  mpz_t two_value;
  mpz_t zero_value;
  mpz_t one_value;
  mpz_init(two_value);
  mpz_init(zero_value);
  mpz_init(one_value);
  mpz_set_str(two_value, "2", 10);
  mpz_set_str(zero_value, "0", 10);
  mpz_set_str(one_value, "1", 10);
  mpz_t modulo_result;
  mpz_init(modulo_result);
  mpz_mod(modulo_result, n, two_value);

  if (mpz_cmp(n, two_value) == 0)
  {
    return true;
  }
  else if (mpz_cmp(modulo_result, zero_value))
  {
    cout << "HHHHHHH " << endl;
    return false;
  }

  mpz_t n1;
  mpz_init(n1);
  mpz_t k;
  mpz_init(k);
  mpz_t u;
  mpz_init(u);
  mpz_set_str(k, "0", 10);
  mpz_sub(n1, n, one_value);
  mpz_set(u, n1);

  mpz_mod(modulo_result, u, two_value);
  /*
  while (mpz_cmp(modulo_result, zero_value))
  {

  }*/





  /*
  cout << "GGGG" << endl;
  mpz_t two_value;
  mpz_set_str(two_value, "2", 10);
  mpz_t modulo_result;
  //mpz_init(modulo_result,NULL);
  cout << mpz_mod(modulo_result, n, two_value) << endl;

  if (mpz_cmp(n, two_value) == 0)
  {
    return true;
  }

  else if (mpz_mod(modulo_result, n, two_value) == true)
  {
    cout << "HHHHHHH " << endl;
    return false;
  }
  //mpz_t n1;
  //mpz_set_ui(n1, n_value - 1);
  //cout << mpz_get_ui(n1) << endl;*/
  return true;




}

int main (int argc, char **argv) {
  int opt;
  bool g = false;
  int modulus_length;
  //gmp_randstate_t state;
  //mpz_t seed;
  unsigned long int seed_value;

  while ((opt = getopt (argc, argv, ":g:")) != -1)
  {
  switch (opt)
  {
    case 'g':
      modulus_length = atoi(optarg);
      cout << modulus_length << endl;
      g = true;
    break;



    default:
      cerr << "Error - Bad parametres!\n";
      exit(EXIT_FAILURE);
    break;

  }
  }

  // Create new random number generator state, and initialize state with the Mersenne Twister algorithm.
  /*mpz_class ran;
  gmp_randclass rr(gmp_randinit_default);
  rr.seed(100000U);
  ran =rr.get_z_bits(125);
  long int random=ran.get_ui();*/

  gmp_randclass r(gmp_randinit_default);
  r.seed(time(NULL));
  mpz_class random_number;
  mpz_t msb_bit;
  mpz_t random_number_value;
  //mpz_randclass ran;
//  gmp_randinit_mt(state);

  // Seed random number generator.
  //mpz_t seed;
  //mpz_init_set_ui(seed, 100000U);
  //gmp_randseed_ui(state, seed);
  random_number = r.get_z_bits(modulus_length);
  long int random = random_number.get_ui();
  cout << random << endl;
  //mpz_init(msb_bit);
//  mpz_set_ui(msb_bit, modulus_length - 1);
  //cout << mpz_get_ui(msb_bit) << endl;

  mpz_init(random_number_value);
  mpz_set_ui(random_number_value, random);
  cout << "FFFF" << endl;
  mpz_setbit(random_number_value, modulus_length - 1);
  cout << "FFFF" << endl;
  cout << mpz_get_ui(random_number_value) << endl;
  //cout << "FFFF" << endl;
  miller_rabin_test(random_number_value);











    return 0;
}
