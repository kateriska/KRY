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
  mpz_set_ui(random_number_value, random);
  mpz_setbit(random_number_value, modulus_length - 1);
  cout << mpz_get_ui(random_number_value) << endl;

  






    return 0;
}
