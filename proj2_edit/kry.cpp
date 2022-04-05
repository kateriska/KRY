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


bool miller_rabin_test(mpz_t n, int modulus_length)
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

  if (mpz_cmp(n, two_value) < 0)
  {
    return false;
  }
  if (mpz_cmp(n, two_value) != 0 && mpz_cmp(modulo_result, zero_value) == 0)
  {
    return false;
  }

  mpz_t k;
  mpz_init(k);
  mpz_sub(k, n, one_value);

  mpz_t modulo_result_k;
  mpz_init(modulo_result_k);

  mpz_mod(modulo_result_k, k, two_value);

  if (mpz_cmp(modulo_result_k, zero_value) == 0)
  {
    while(true)
    {
      mpz_fdiv_q(k,k,two_value);

      mpz_t modulo_result_k_iteration;
      mpz_init(modulo_result_k_iteration);

      mpz_mod(modulo_result_k_iteration, k, two_value);

      if (mpz_cmp(modulo_result_k_iteration, zero_value) == 0)
      {
        continue;
      }
      else
      {
        break;
      }
    }
  }
  cout << "Counted k" << endl;
  mpz_out_str(stdout,10,k);





  for (int i = 0; i < 5; i++)
  {
    //unsigned long int a_int = rand() % generated_a_range + 2;
    //cout << a_int << endl;
    mpz_t a;
    mpz_init(a);
    gmp_randstate_t state;
    gmp_randinit_mt (state);
    gmp_randseed_ui(state, time(NULL));
    mpz_urandomm(a,state,n);
    cout << "Getting random number a" <<  endl;
    //mpz_set_str(a, "39", 10);
    mpz_out_str(stdout,10,a);

    //cout << "DDDD" <<  endl;
    mpz_t k_copy;
    mpz_init(k_copy);
    mpz_set(k_copy, k);

    //mpz_out_str(stdout,10,k_copy);

    // counting modular exponentiation
    // (modular_exponentiation_base^modular_exponentiation_exponent) % modular_exponentiation_modulo
    mpz_t modular_exponentiation_base;
    mpz_t modular_exponentiation_exponent;
    mpz_t modular_exponentiation_modulo;

    mpz_init(modular_exponentiation_base);
    mpz_init(modular_exponentiation_exponent);
    mpz_init(modular_exponentiation_modulo);

    mpz_set(modular_exponentiation_base, a);
    //cout << "Getting modular_exponent_base" << endl;
    //mpz_out_str(stdout,10,modular_exponentiation_base);
    mpz_set(modular_exponentiation_exponent, k_copy);
    mpz_set(modular_exponentiation_modulo, n);

    mpz_t x;
    mpz_init(x);
    mpz_set_str(x, "1", 10);
    mpz_t y;
    mpz_init(y);
    mpz_set(y, modular_exponentiation_base);
    //cout << "ZZZ" << endl;




      while (mpz_cmp(modular_exponentiation_exponent, zero_value) > 0)
      {
        mpz_t modulo_exponent;
        mpz_init(modulo_exponent);
        mpz_mod(modulo_exponent, modular_exponentiation_exponent, two_value);
        if (mpz_cmp(modulo_exponent, one_value) == 0)
        {
          mpz_t mul_operation;
          mpz_init(mul_operation);
          mpz_mul(mul_operation, x, y);
          //mpz_out_str(stdout,10,y);
          //cout << "Update x" << endl;
          mpz_mod(x, mul_operation, modular_exponentiation_modulo);
        }
        mpz_t mul_operation;
        mpz_init(mul_operation);
        mpz_mul(mul_operation, y, y);
        mpz_mod(y, mul_operation, modular_exponentiation_modulo);
        mpz_fdiv_q(modular_exponentiation_exponent, modular_exponentiation_exponent, two_value);
        //mpz_out_str(stdout,10,exponent);





    }
    mpz_t result;
    mpz_init(result);
    mpz_mod(result, x, modular_exponentiation_modulo);

    cout << "Result of modular exponentiation" <<  endl;
    mpz_out_str(stdout,10,result);

    // end of counting modular exponentiation
    mpz_t n1;
    mpz_init(n1);
    mpz_sub(n1, n, one_value);

    // couning of (multiply_modulo_factor1 * multiply_modulo_factor2) % multiply_modulo_modulo
    while (mpz_cmp(k_copy, n1) != 0 && mpz_cmp(result, one_value) != 0 && mpz_cmp(result, n1) != 0)
    {
      mpz_t multiply_modulo_factor1;
      mpz_t multiply_modulo_factor2;
      mpz_t multiply_modulo_modulo;
      mpz_init(multiply_modulo_factor1);
      mpz_init(multiply_modulo_factor2);
      mpz_init(multiply_modulo_modulo);

      mpz_set(multiply_modulo_factor1, result);
      mpz_set(multiply_modulo_factor2, result);
      mpz_set(multiply_modulo_modulo, n);

      mpz_t x;
      mpz_init(x);
      mpz_set_str(x, "0", 10);
      mpz_t y;
      mpz_init(y);
      mpz_mod(y, multiply_modulo_factor1, multiply_modulo_modulo);

      while (mpz_cmp(multiply_modulo_factor2, zero_value) > 0)
      {
        mpz_t modulo_factor2;
        mpz_init(modulo_factor2);
        mpz_mod(modulo_factor2, multiply_modulo_factor2, two_value);

        if (mpz_cmp(modulo_factor2, one_value) == 0)
        {
          mpz_t add_operation;
          mpz_init(add_operation);
          mpz_add(add_operation, x, y);
          mpz_mod(x, add_operation, multiply_modulo_modulo);
        }

        mpz_t mul_operation;
        mpz_init(mul_operation);
        mpz_mul(mul_operation, y, two_value);
        mpz_mod(y, mul_operation, multiply_modulo_modulo);
        mpz_fdiv_q(multiply_modulo_factor2, multiply_modulo_factor2, two_value);

      }

      mpz_mod(result, x, multiply_modulo_modulo);
        // end of multiply modulo operation
      mpz_mul(k_copy, k_copy, two_value);

    }



    cout << "Result of multiply_modulo" <<  endl;
    mpz_out_str(stdout,10,result);
    //cout << "Res2" <<  endl;

    mpz_t k_copy_modulo;
    mpz_init(k_copy_modulo);
    mpz_mod(k_copy_modulo, k_copy, two_value);
    if (mpz_cmp(result, n1) != 0 && mpz_cmp(k_copy_modulo, zero_value) == 0)
    {
      return false;
    }



  }

  cout << "Potential prime" << endl;
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



  gmp_randclass r(gmp_randinit_default);
  r.seed(time(NULL));
  mpz_class random_number;
  mpz_t msb_bit;
  mpz_t random_number_value;
  mpz_t modulus_length_value;

  random_number = r.get_z_bits(modulus_length);
  mpz_init(random_number_value);
  mpz_init(modulus_length_value);

  mpz_set(random_number_value, random_number.get_mpz_t());
  size_t size = mpz_sizeinbase(random_number_value, 2);
  cout << size << endl;
  mpz_init(msb_bit);
  mpz_set_ui(msb_bit, modulus_length - 1);
  mpz_out_str(stdout,10,random_number_value);
  mpz_set_ui(modulus_length_value, modulus_length);

  mpz_t prime_value;
  mpz_init(prime_value);
  mpz_set_str(prime_value, "97", 10);

  mpz_t one_value;
  mpz_init(one_value);
  mpz_set_str(one_value, "1", 10);

  //miller_rabin_test(prime_value, modulus_length);


  vector <unsigned long int> potential_primes;
  int generated_prime_count = 2;
  while (generated_prime_count > 0)
  {
  if (miller_rabin_test(random_number_value, modulus_length))
  {
    //cout << mpz_get_ui(random_number_value) << " is potential prime" << endl;
    cout << "This number is potential prime:" << endl;
    mpz_out_str(stdout,10,random_number_value);
    generated_prime_count--;
    unsigned long int random_number_value_int = mpz_get_ui(random_number_value);
    potential_primes.push_back(random_number_value_int);

  }
  mpz_add(random_number_value, random_number_value, one_value);
  }

  cout << "Prime numms:" << endl;
  cout << potential_primes.at(0) << endl;
  cout << potential_primes.at(1) << endl;




    return 0;
}
