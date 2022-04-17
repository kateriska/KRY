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

// test for primality of potential prime candidate
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

  // even numbers except 2 can't be primes
  if (mpz_cmp(n, two_value) < 0)
  {
    mpz_clear(two_value);
    mpz_clear(zero_value);
    mpz_clear(one_value);
    mpz_clear(modulo_result);
    return false;
  }
  if (mpz_cmp(n, two_value) != 0 && mpz_cmp(modulo_result, zero_value) == 0)
  {
    mpz_clear(two_value);
    mpz_clear(zero_value);
    mpz_clear(one_value);
    mpz_clear(modulo_result);
    return false;
  }

  mpz_t k;
  mpz_init(k);
  mpz_sub(k, n, one_value);

  mpz_t modulo_result_k;
  mpz_init(modulo_result_k);

  mpz_mod(modulo_result_k, k, two_value);

  mpz_t modulo_result_k_iteration;
  mpz_init(modulo_result_k_iteration);

  if (mpz_cmp(modulo_result_k, zero_value) == 0)
  {
    while(true)
    {
      mpz_fdiv_q(k,k,two_value);

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

  // declaration of variables outside a loop of testing witnesses for simple cleaning after processing
  mpz_t a; // randomly generated witness
  mpz_init(a);
  mpz_t k_copy;
  mpz_init(k_copy);

  // variables for modular exponentiation
  mpz_t modular_exponentiation_base;
  mpz_t modular_exponentiation_exponent;
  mpz_t modular_exponentiation_modulo;

  mpz_init(modular_exponentiation_base);
  mpz_init(modular_exponentiation_exponent);
  mpz_init(modular_exponentiation_modulo);

  mpz_t x;
  mpz_init(x);
  mpz_t y;
  mpz_init(y);

  mpz_t modulo_exponent;
  mpz_init(modulo_exponent);
  mpz_t mul_operation;
  mpz_init(mul_operation);

  mpz_t result;
  mpz_init(result);

  mpz_t n1;
  mpz_init(n1);

  // variables for multiply and modulo operation
  mpz_t multiply_modulo_factor1;
  mpz_t multiply_modulo_factor2;
  mpz_t multiply_modulo_modulo;
  mpz_init(multiply_modulo_factor1);
  mpz_init(multiply_modulo_factor2);
  mpz_init(multiply_modulo_modulo);

  mpz_t modulo_factor2;
  mpz_init(modulo_factor2);

  mpz_t add_operation;
  mpz_init(add_operation);

  mpz_t k_copy_modulo;
  mpz_init(k_copy_modulo);

  gmp_randstate_t state; // random state init for generating possible witnesses
  gmp_randinit_mt (state);

  // find possible witnesses
  for (int i = 0; i < 15; i++)
  {
    mpz_sub(n1, n, one_value);
    gmp_randseed_ui(state, time(NULL));
    mpz_urandomm(a,state,n1); // generate random number less or equal to n - 2

    mpz_set(k_copy, k);

    // counting modular exponentiation
    // (modular_exponentiation_base^modular_exponentiation_exponent) % modular_exponentiation_modulo

    mpz_set(modular_exponentiation_base, a);
    mpz_set(modular_exponentiation_exponent, k_copy);
    mpz_set(modular_exponentiation_modulo, n);

    mpz_set_str(x, "1", 10);
    mpz_set(y, modular_exponentiation_base);

      while (mpz_cmp(modular_exponentiation_exponent, zero_value) > 0)
      {
        mpz_mod(modulo_exponent, modular_exponentiation_exponent, two_value);
        if (mpz_cmp(modulo_exponent, one_value) == 0)
        {
          mpz_mul(mul_operation, x, y);
          mpz_mod(x, mul_operation, modular_exponentiation_modulo);
        }

        mpz_mul(mul_operation, y, y);
        mpz_mod(y, mul_operation, modular_exponentiation_modulo);
        mpz_fdiv_q(modular_exponentiation_exponent, modular_exponentiation_exponent, two_value);
    }

    mpz_mod(result, x, modular_exponentiation_modulo);
    // end of counting modular exponentiation

    // counting of (multiply_modulo_factor1 * multiply_modulo_factor2) % multiply_modulo_modulo
    while (mpz_cmp(k_copy, n1) != 0 && mpz_cmp(result, one_value) != 0 && mpz_cmp(result, n1) != 0)
    {
      mpz_set(multiply_modulo_factor1, result);
      mpz_set(multiply_modulo_factor2, result);
      mpz_set(multiply_modulo_modulo, n);

      mpz_set_str(x, "0", 10);
      mpz_mod(y, multiply_modulo_factor1, multiply_modulo_modulo);

      while (mpz_cmp(multiply_modulo_factor2, zero_value) > 0)
      {
        mpz_mod(modulo_factor2, multiply_modulo_factor2, two_value);

        if (mpz_cmp(modulo_factor2, one_value) == 0)
        {
          mpz_add(add_operation, x, y);
          mpz_mod(x, add_operation, multiply_modulo_modulo);
        }

        mpz_mul(mul_operation, y, two_value);
        mpz_mod(y, mul_operation, multiply_modulo_modulo);
        mpz_fdiv_q(multiply_modulo_factor2, multiply_modulo_factor2, two_value);

      }

      mpz_mod(result, x, multiply_modulo_modulo);
        // end of multiply modulo operation
      mpz_mul(k_copy, k_copy, two_value);

    }


    mpz_mod(k_copy_modulo, k_copy, two_value);
    if (mpz_cmp(result, n1) != 0 && mpz_cmp(k_copy_modulo, zero_value) == 0)
    {
      // necessary cleaning to avoid memory leaks
      mpz_clear(two_value);
      mpz_clear(zero_value);
      mpz_clear(one_value);

      mpz_clear(modulo_result);

      mpz_clear(k);
      mpz_clear(modulo_result_k);
      mpz_clear(modulo_result_k_iteration);

      mpz_clear(a);
      mpz_clear(k_copy);

      mpz_clear(modular_exponentiation_base);
      mpz_clear(modular_exponentiation_exponent);
      mpz_clear(modular_exponentiation_modulo);

      mpz_clear(x);
      mpz_clear(y);

      mpz_clear(modulo_exponent);
      mpz_clear(mul_operation);

      mpz_clear(result);

      mpz_clear(n1);

      mpz_clear(multiply_modulo_factor1);
      mpz_clear(multiply_modulo_factor2);
      mpz_clear(multiply_modulo_modulo);

      mpz_clear(modulo_factor2);

      mpz_clear(add_operation);
      mpz_clear(k_copy_modulo);

      gmp_randclear(state);
      return false;
    }



  }

  // necessary cleaning to avoid memory leaks
  mpz_clear(two_value);
  mpz_clear(zero_value);
  mpz_clear(one_value);

  mpz_clear(modulo_result);

  mpz_clear(k);
  mpz_clear(modulo_result_k);
  mpz_clear(modulo_result_k_iteration);

  mpz_clear(a);
  mpz_clear(k_copy);

  mpz_clear(modular_exponentiation_base);
  mpz_clear(modular_exponentiation_exponent);
  mpz_clear(modular_exponentiation_modulo);

  mpz_clear(x);
  mpz_clear(y);

  mpz_clear(modulo_exponent);
  mpz_clear(mul_operation);

  mpz_clear(result);

  mpz_clear(n1);

  mpz_clear(multiply_modulo_factor1);
  mpz_clear(multiply_modulo_factor2);
  mpz_clear(multiply_modulo_modulo);

  mpz_clear(modulo_factor2);

  mpz_clear(add_operation);
  mpz_clear(k_copy_modulo);

  gmp_randclear(state);

  return true;
}


void compute_gcd(mpz_t out,mpz_t a, mpz_t b)
{
  mpz_t zero_value;
  mpz_init(zero_value);
  mpz_set_str(zero_value, "0", 10);

  mpz_t modulo_operation;
  mpz_init(modulo_operation);

  mpz_mod(modulo_operation, a, b);

  if (mpz_cmp(modulo_operation, zero_value) == 0)
  {
      mpz_set(out,b);
      mpz_clear(modulo_operation);
      mpz_clear(zero_value);
      return;
  }

  compute_gcd(out, b, modulo_operation);
  mpz_clear(modulo_operation);
  mpz_clear(zero_value);
  return;
}

/*
Veřejný modulus nejdříve zkontrolujte pomocí metody '''triviálního (zkusmého) dělení''' pro '''prvních 1 000 000 čísel'''.
Metoda triviálního dělení je nejjednodušší metodou pro faktorizaci celých čísel
Metoda pracuje na základě zvolení počátečního dělitele (např. 2) a následném ověření, zda dělitel opravdu dělí zadané číslo beze zbytku.
Pokud ne, dojde k inkrementaci hodnoty dělitele a opět se ověřuje dělitelnost beze zbytku.
Tento proces je opakován tak dlouho, dokud není nalezen správný dělitel.
*/
void trivial_division(mpz_t out, mpz_t n)
{
  mpz_t divisor;
  mpz_init(divisor);
  mpz_set_str(divisor, "2", 10);

  mpz_t divisor_end;
  mpz_init(divisor_end);
  mpz_set_str(divisor_end, "1000002", 10);

  mpz_t zero_value;
  mpz_init(zero_value);
  mpz_set_str(zero_value, "0", 10);

  mpz_t one_value;
  mpz_init(one_value);
  mpz_set_str(one_value, "1", 10);

  mpz_t modulo_result;
  mpz_init(modulo_result);

  mpz_mod(modulo_result, n, divisor);

  mpz_t modulo_result_iteration;
  mpz_init(modulo_result_iteration);
  if (mpz_cmp(modulo_result, zero_value) != 0)
  {
    while (true)
    {
      mpz_mod(modulo_result_iteration, n, divisor);

      // found correct divisor
      if (mpz_cmp(modulo_result_iteration, zero_value) == 0)
      {
        mpz_set(out, divisor);

        mpz_clear(divisor);
        mpz_clear(divisor_end);
        mpz_clear(zero_value);
        mpz_clear(one_value);
        mpz_clear(modulo_result);
        mpz_clear(modulo_result_iteration);
        return;
      }

      mpz_add(divisor, divisor, one_value);

      // only test first 1 000 000 divisors, so if finding divisor was stil unsucessful, set divisor to 1 and possibly find them with Fermat factorization
      if (mpz_cmp(divisor, divisor_end) == 0)
      {
        mpz_set(out, one_value);

        mpz_clear(divisor);
        mpz_clear(divisor_end);
        mpz_clear(zero_value);
        mpz_clear(one_value);
        mpz_clear(modulo_result);
        mpz_clear(modulo_result_iteration);
        return;
      }
    }
  }

  mpz_set(out, divisor);
  mpz_clear(divisor);
  mpz_clear(divisor_end);
  mpz_clear(zero_value);
  mpz_clear(one_value);
  mpz_clear(modulo_result);
  mpz_clear(modulo_result_iteration);
  return;
}

void fermat_factorization(mpz_t out, mpz_t n)
{
  mpz_t a;
  mpz_init(a);

  mpz_t n_sqrt;
  mpz_init(n_sqrt);

  mpz_t one_value;
  mpz_init(one_value);
  mpz_set_str(one_value, "1", 10);

  mpz_sqrt(n_sqrt, n);
  mpz_add(a, n_sqrt, one_value);

  mpz_t mul_operation;
  mpz_init(mul_operation);
  mpz_mul(mul_operation, a, a);

  if (mpz_cmp(mul_operation, n) == 0)
  {
    mpz_set(out, a);

    mpz_clear(a);
    mpz_clear(n_sqrt);
    mpz_clear(one_value);
    mpz_clear(mul_operation);
    return;
  }

  mpz_t b;
  mpz_init(b);

  mpz_t b2;
  mpz_init(b2);

  mpz_t mul_operation_iteration_a;
  mpz_init(mul_operation_iteration_a);
  mpz_t mul_operation_iteration_b;
  mpz_init(mul_operation_iteration_b);

  while (true)
  {
    mpz_mul(mul_operation_iteration_a, a, a);
    mpz_sub(b2, mul_operation_iteration_a, n);
    mpz_sqrt(b,b2);

    mpz_mul(mul_operation_iteration_b, b, b);
    if (mpz_cmp(mul_operation_iteration_b, b2) == 0)
    {
      break;
    }
    else
    {
      mpz_add(a, a, one_value);
    }

  }

  mpz_t p;
  mpz_init(p);
  mpz_t q;
  mpz_init(q);

  mpz_t computed_n;
  mpz_init(computed_n);

  mpz_sub(p, a, b);
  mpz_add(q, a, b);
  mpz_set(out, p);
  mpz_mul(computed_n, p, q);

  /*
  // check of correctness
  if (mpz_cmp(n, computed_n) == 0)
  {
    cout << "Fermat factorization was success" << endl;
  }
  */

  mpz_clear(a);
  mpz_clear(n_sqrt);
  mpz_clear(one_value);
  mpz_clear(mul_operation);

  mpz_clear(b);
  mpz_clear(b2);
  mpz_clear(mul_operation_iteration_a);
  mpz_clear(mul_operation_iteration_b);

  mpz_clear(p);
  mpz_clear(q);

  mpz_clear(computed_n);

  return;
}



int main (int argc, char **argv) {
  int opt;
  bool g = false;
  bool e = false;
  bool d = false;
  bool b = false;
  int modulus_length;

  while ((opt = getopt (argc, argv, ":g:e:d:b:")) != -1)
  {
  switch (opt)
  {
    case 'g':
      modulus_length = atoi(optarg);
      g = true;
    break;

    case 'e':
      e = true;
    break;

    case 'd':
      d = true;
    break;

    case 'b':
      b = true;
    break;

    default:
      cerr << "Error - Bad parametres!\n";
      exit(EXIT_FAILURE);
    break;

  }
  }


  if (g == true) {
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

  mpz_init(msb_bit);
  // set msb bit to 1
  mpz_set_ui(msb_bit, modulus_length - 1);
  mpz_set_ui(modulus_length_value, modulus_length);

  mpz_t one_value;
  mpz_init(one_value);
  mpz_set_str(one_value, "1", 10);

  mpz_t zero_value;
  mpz_init(zero_value);
  mpz_set_str(zero_value, "0", 10);

  int generated_prime_count = 2;

  mpz_t p;
  mpz_init(p);
  mpz_t q;
  mpz_init(q);

  mpz_t p1;
  mpz_init(p1);
  mpz_t q1;
  mpz_init(q1);

  // generate 2 potential primes
  while (generated_prime_count > 0)
  {
  if (miller_rabin_test(random_number_value, modulus_length))
  {
    if (generated_prime_count == 2)
    {
      mpz_set(p, random_number_value);
    }
    else if (generated_prime_count == 1)
    {
      mpz_set(q, random_number_value);
    }

    generated_prime_count--;

  }
  mpz_add(random_number_value, random_number_value, one_value);
  }

  mpz_t n;
  mpz_init(n);

  mpz_t phi_n;
  mpz_init(phi_n);

  mpz_mul(n, p, q);

  mpz_sub(p1, p, one_value);
  mpz_sub(q1, q, one_value);

  mpz_mul(phi_n, p1, q1);

  // Zvol náhodně e mezi 1 a phi( n ) tak, že gcd(e, phi( n )) = 1
  mpz_t e;
  mpz_init(e);

  gmp_randstate_t state;
  gmp_randinit_mt (state);
  gmp_randseed_ui(state, time(NULL));
  mpz_urandomm(e,state,phi_n);

  mpz_t gcd_res;
  mpz_init(gcd_res);
  while (true)
  {
    compute_gcd(gcd_res, e, phi_n);

    if (mpz_cmp(gcd_res, one_value) == 0)
    {
      break;
    }
    mpz_add(e, e, one_value);
  }

  // compute of multiplicative inverse
  mpz_t a;
  mpz_init(a);
  mpz_t m;
  mpz_init(m);

  mpz_set(a, e);
  mpz_set(m, phi_n);

  mpz_t m0;
  mpz_init(m0);
  mpz_set(m0, m);
  mpz_t y;
  mpz_init(y);
  mpz_t x;
  mpz_init(x);
  mpz_set_str(y, "0", 10);
  mpz_set_str(x, "1", 10);

  mpz_t quotient;
  mpz_init(quotient);
  mpz_t t;
  mpz_init(t);
  mpz_t mul_operation;
  mpz_init(mul_operation);


  if (mpz_cmp(m, one_value) == 0)
  {
    mpz_set(x, zero_value);
  }
  else {
  while (mpz_cmp(a, one_value) > 0)
  {
    mpz_fdiv_q(quotient, a, m);

    mpz_set(t,m);

    mpz_mod(m, a, m);
    mpz_set(a,t);
    mpz_set(t,y);

    mpz_mul(mul_operation, quotient, y);
    mpz_sub(y, x, mul_operation);
    mpz_set(x,t);

  }
}

  if (mpz_cmp(x, zero_value) < 0)
  {
    mpz_add(x,x,m0);
  }

  cout << modulus_length << " ";
  cout << "0x";
  mpz_out_str(stdout,16,p);
  cout << " ";
  cout << "0x";
  mpz_out_str(stdout,16,q);
  cout << " ";
  cout << "0x";
  mpz_out_str(stdout,16,n);
  cout << " ";
  cout << "0x";
  mpz_out_str(stdout,16,e);
  cout << " ";
  cout << "0x";
  mpz_out_str(stdout,16,x);

  putchar('\n');

  /*
  cout << "Modulus length: " << endl;
  cout << modulus_length;
  putchar('\n');
  cout << "P: " << endl;
  mpz_out_str(stdout,10,p);
  putchar('\n');
  cout << "Q: " << endl;
  mpz_out_str(stdout,10,q);
  putchar('\n');
  cout << "N: " << endl;
  mpz_out_str(stdout,10,n);
  putchar('\n');
  cout << "E: " << endl;
  mpz_out_str(stdout,10,e);
  putchar('\n');
  cout << "X: " << endl;
  mpz_out_str(stdout,10,x);
  putchar('\n');
  */
  // cleaning to avoid memory leaks
  mpz_clear(random_number_value);
  mpz_clear(modulus_length_value);

  mpz_clear(msb_bit);

  mpz_clear(one_value);
  mpz_clear(zero_value);

  mpz_clear(p);
  mpz_clear(q);

  mpz_clear(p1);
  mpz_clear(q1);

  mpz_clear(n);
  mpz_clear(phi_n);
  mpz_clear(e);

  mpz_clear(gcd_res);

  mpz_clear(a);
  mpz_clear(m);
  mpz_clear(m0);
  mpz_clear(y);
  mpz_clear(x);
  mpz_clear(quotient);
  mpz_clear(t);
  mpz_clear(mul_operation);

  gmp_randclear(state);



}

/*
'''Šifrování''' (0.4b)

vstup: ./kry -e E N M

výstup: C

'''Dešifrování''' (0.4b)

vstup: ./kry -d D N C
*/

if (e == true)
{
  mpz_t e;
  mpz_init(e);
  mpz_t n;
  mpz_init(n);
  mpz_t m;
  mpz_init(m);

  mpz_set_str(e, argv[2], 0);
  mpz_set_str(n, argv[3], 0);
  mpz_set_str(m, argv[4], 0);

  mpz_t c;
  mpz_init(c);

  mpz_powm(c, m, e, n);

  cout << "0x";
  mpz_out_str(stdout,16,c);
  putchar('\n');

  mpz_clear(e);
  mpz_clear(n);
  mpz_clear(m);
  mpz_clear(c);
}

if (d == true)
{
  mpz_t d;
  mpz_init(d);
  mpz_t n;
  mpz_init(n);
  mpz_t c;
  mpz_init(c);

  mpz_set_str(d, argv[2], 0);
  mpz_set_str(n, argv[3], 0);
  mpz_set_str(c, argv[4], 0);

  mpz_t m;
  mpz_init(m);

  mpz_powm(m, c, d, n);

  cout << "0x";
  mpz_out_str(stdout,16,m);
  putchar('\n');

  mpz_clear(d);
  mpz_clear(n);
  mpz_clear(c);
  mpz_clear(m);
}

/*
'''Prolomení RSA''' (2.6b)

vstup: ./kry -b N

výstup: P
*/

if (b == true)
{
  mpz_t one_value;
  mpz_init(one_value);
  mpz_set_str(one_value, "1", 10);

  mpz_t n;
  mpz_init(n);
  mpz_set_str(n, argv[2], 0);

  mpz_t trivial_division_res;
  mpz_init(trivial_division_res);
  trivial_division(trivial_division_res, n);

  if (mpz_cmp(trivial_division_res, one_value) != 0)
  {
    cout << "0x";
    mpz_out_str(stdout,16,trivial_division_res);
    putchar('\n');
  }

  mpz_t fermat_factorization_res;
  mpz_init(fermat_factorization_res);

  // trivial division unsucessful, need to call more effective fermat factorization
  if (mpz_cmp(trivial_division_res, one_value) == 0)
  {
    //cout << "Need to do fermat factorization" << endl;
    fermat_factorization(fermat_factorization_res, n);
    cout << "0x";
    mpz_out_str(stdout,16,fermat_factorization_res);
    putchar('\n');
  }

  mpz_clear(one_value);
  mpz_clear(n);
  mpz_clear(trivial_division_res);
  mpz_clear(fermat_factorization_res);

}

  return 0;
}
