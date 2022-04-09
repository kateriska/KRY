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

unsigned long int modular_exponent(mpz_t a, mpz_t u, mpz_t n);
bool find_witness(mpz_t a, mpz_t u, mpz_t n, mpz_t n1, mpz_t k);

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

  else if (mpz_cmp(modulo_result, zero_value) == 0)
  {
    cout << "Odd number - cant be prime" << endl;
    return false;
  }


  mpz_t n1;
  mpz_init(n1); // p1
  mpz_t k;
  mpz_init(k);
  mpz_t u;
  mpz_init(u);
  mpz_set_str(k, "0", 10); // u
  mpz_sub(n1, n, one_value);
  mpz_set(u, n1); // r = u
  cout << "Printing u " << endl;
  cout << mpz_get_ui(u) << endl;
  mpz_t modulo_result_shift;
  mpz_init(modulo_result_shift);

  mpz_mod(modulo_result_shift, u, two_value);

  if (mpz_cmp(modulo_result_shift, zero_value) == 0) {
  while (true)
  {
    cout << "DKDKD" << endl;
    mpz_fdiv_q_2exp(u,u,1); // shift right'
    //unsigned long int u_int = mpz_get_ui(u);
    //u_int = u_int >> 1;
    //mpz_set_ui(u, u_int);

    //cout << "Printing u" << endl;
    //cout << u_int << endl;
    mpz_add(k, k, one_value);
    mpz_t modulo_result_shift_iteration;
    mpz_init(modulo_result_shift_iteration);

    mpz_mod(modulo_result_shift_iteration, u, two_value);
    if (mpz_cmp(modulo_result_shift_iteration, zero_value) == 0)
    {
      continue;
    }
    else
    {
      break;
    }
  }
}
 // ok nahoru
  /*
  mpz_t left_side;
  mpz_t right_side;
  mpz_init(left_side);
  mpz_init(right_side);
  mpz_sub(left_side, n, one_value);

  unsigned long int k_int = mpz_get_ui(k);
  cout << "Printing k int" << endl;
  cout << k_int << endl;

  mpz_t right_side_power;
  mpz_init(right_side_power);
  mpz_pow_ui(right_side_power, two_value, k_int);
  mpz_mul(right_side, right_side_power, u);

  unsigned long int left_side_int = mpz_get_ui(left_side);
  unsigned long int right_side_int = mpz_get_ui(right_side);

  cout << "Left and right side " << endl;
  cout << right_side_int << endl;
  cout << left_side_int << endl;
  */


  mpz_t num1;
  mpz_init(num1);
  mpz_t num2;
  mpz_init(num2);
  mpz_t num3;
  mpz_init(num3);

  mpz_set_str(num1, "3", 10);
  mpz_set_str(num2, "2", 10);
  mpz_set_str(num3, "2", 10);

  //cout << modular_exponent(num1, num2, num3) << endl;

  //unsigned long int generated_a_range = mpz_get_ui(n) - 2 + 1;
  mpz_t generated_a_range;
  mpz_init(generated_a_range);

  mpz_sub(generated_a_range, n, two_value);
  mpz_add(generated_a_range, generated_a_range, 1);
  cout << "LLLL" << endl;
  cout << generated_a_range << endl;

  mpz_t a,temp;
  mpz_init(a);
  mpz_init(temp);
  gmp_randstate_t state;
  gmp_randinit_mt (state);
  mpz_urandomb (temp, state,lambda);
  mpz_nextprime (group_size, temp);

  //srand(time(0));
  //unsigned long int a_int = rand() % generated_a_range + 2;
  //cout << a_int << endl;

  // ok nahoru

  //mpz_t a;
  //mpz_init(a);
  //mpz_set_ui(a, a_int);






  /*
  for (int i = 0; i < 5; i++)
  {
    unsigned long int a_int = rand() % generated_a_range + 2;
    cout << a_int << endl;
    mpz_t a;
    mpz_init(a);
    mpz_set_ui(a, a_int);
    if (find_witness(a, u, n, n1, k))
    {
      cout << "Witness that p is not prime " << endl;
      return false;
    }
  }
  */

  return true;
}

bool find_witness(mpz_t a, mpz_t u, mpz_t n, mpz_t n1, mpz_t k)
{
  mpz_t two_value;
  mpz_init(two_value);
  mpz_set_str(two_value, "2", 10);

  mpz_t three_value;
  mpz_init(three_value);
  mpz_set_str(three_value, "3", 10);

  cout << "===PARAMS===" << endl;
  cout << mpz_get_ui(a) << endl;
  cout << mpz_get_ui(u) << endl;
  cout << mpz_get_ui(n) << endl;
  unsigned long int modular_exponent_result = modular_exponent(a, u, n);
  cout << "Modular exponent result " << modular_exponent_result << endl;


  unsigned long int n_int = mpz_get_ui(n);
  if (modular_exponent_result == 1)
  {
    return false;
  }

  unsigned long int k_int = mpz_get_ui(k);
  cout << "K " << k_int << endl;

  for (unsigned long int i = 0; i < k_int; i++)
  {
    cout << "I " << i << endl;
  //  modular_exponent_result = modular_exponent(a, u, n);

    mpz_t modular_exponent_multiply;
    mpz_t modular_exponent_pow;
    mpz_init(modular_exponent_multiply);
    mpz_init(modular_exponent_pow);
    mpz_pow_ui(modular_exponent_pow, two_value, i);
    cout << "Modular exponent power" << endl;
    cout << mpz_get_ui(modular_exponent_pow) << endl;
    mpz_mul(modular_exponent_multiply, modular_exponent_pow, u);
    cout << "Modular exponent multiply" << endl;
    cout << mpz_get_ui(modular_exponent_multiply) << endl;
    unsigned long int modular_exponent_result_iteration = modular_exponent(a, modular_exponent_multiply, n);

    cout << "N1 " <<  mpz_get_ui(n1) << endl;
    cout << "New modular exponent result " << modular_exponent_result_iteration << endl;
    if (modular_exponent_result_iteration == mpz_get_ui(n1))
    {
      return false;
    }

  }


  return true;
}

// https://www.geeksforgeeks.org/count-total-bits-number/
unsigned int countBits(unsigned int n)
{
   unsigned int count = 0;
   while (n)
   {
        count++;
        n >>= 1;
    }
    return count;
}
// https://www.geeksforgeeks.org/binary-representation-of-a-given-number/
/*
string bin(unsigned int n)
{
    if (n > 1)
        bin(n >> 1);

    printf("%d", n & 1);
    int res = n & 1;
    cout << "Counting res " << endl;
    cout << res << endl;
    cout << "Counting res endl" << endl;
    string res_str = to_string(res);
    return res_str;
}
*/
// https://stackoverflow.com/questions/22746429/c-decimal-to-binary-converting
string toBinary(unsigned long int n) {
    if (n==0) return "0";
    else if (n==1) return "1";
    else if (n%2 == 0) return toBinary(n/2) + "0";
    else if (n%2 != 0) return toBinary(n/2) + "1";
}



unsigned long int modular_exponent(mpz_t a, mpz_t u, mpz_t n)
{
  mpz_t two_value;
  mpz_init(two_value);
  mpz_set_str(two_value, "2", 10);

  mpz_t zero_value;
  mpz_init(zero_value);
  mpz_set_str(zero_value, "0", 10);

  mpz_t result;
  mpz_init(result);
  mpz_set_str(result, "1", 10);

  unsigned long int u_int = mpz_get_ui(u);


  cout << "Var " << countBits(u_int) <<  endl;

  string u_bin = toBinary(u_int);
  cout << "VAR n" << u_bin.length() << endl;

  for (int i = 0; i < u_bin.length(); i++)
  {
    cout << "HHHHHHHHHHHH" << endl;
    cout << u_bin[i] << endl;
    cout << "FFFFFFFFFFFF" << endl;

    mpz_pow_ui(result, result, 2);
    if (u_bin[i] == '1')
    {
      mpz_mul(result, result, a);
    }
    mpz_mod(result, result, n);


  }

  //mpz_mod(a, a, n);

  /*
  if (mpz_cmp(u, zero_value) > 0) {
  while (true)
  {
    mpz_t modulo_result_iteration;
    mpz_init(modulo_result_iteration);
    mpz_mod(modulo_result_iteration, u, two_value);
    if (mpz_cmp(modulo_result_iteration, zero_value) != 0) // if u is odd
    {
      mpz_t odd_multiply;
      mpz_t odd_modulo;
      mpz_init(odd_multiply);
      mpz_init(odd_modulo);
      mpz_mul(odd_multiply, result, a);
      mpz_mod(odd_modulo, odd_multiply, n);
      mpz_set(result, odd_modulo);
    }

    mpz_fdiv_q_2exp(u,u,1); // shift right
    mpz_t a_multiply;
    mpz_init(a_multiply);
    mpz_mul(a_multiply, a, a);
    mpz_mod(a, a_multiply, n);

    if (mpz_cmp(u, zero_value) > 0)
    {
      continue;
    }
    else
    {
      break;
    }
  }
}
*/
return mpz_get_ui(result);





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

  random_number = r.get_z_bits(1024);
  long int random = random_number.get_ui();
  //mpz_out_str(stdout,16,p);
  cout << random << endl;
  mpz_init(random_number_value);
  //mpz_set(random_number_value, random_number);

  //printf("%zu",size);
  mpz_set(random_number_value, random_number.get_mpz_t());
  size_t size = mpz_sizeinbase(random_number_value, 2);
  cout << size << endl;
  mpz_init(msb_bit);
  mpz_set_ui(msb_bit, modulus_length - 1);
  mpz_out_str(stdout,10,random_number_value);


  vector <unsigned long int> potential_primes;
  mpz_t p;
  mpz_init(p);
  mpz_t q;
  mpz_init(q);

  while (generated_prime_count > 0)
  {
  if (miller_rabin_test(random_number_value))
  {
    cout << mpz_get_ui(random_number_value) << " is potential prime" << endl;

    if (generated_prime_count == 2)
    {
      mpz_set(p, random_number_value);
    }
    else if (generated_prime_count == 1)
    {
      mpz_set(q, random_number_value);
    }
    generated_prime_count--;
    unsigned long int random_number_value_int = mpz_get_ui(random_number_value);
    potential_primes.push_back(random_number_value_int);

  }
  mpz_add(random_number_value, random_number_value, one_value);
  }

  cout << "Prime numms:" << endl;
  cout << potential_primes.at(0) << endl;
  cout << potential_primes.at(1) << endl;

  mpz_t n;
  mpz_init(n);
  mpz_mul(n, p, q);

















    return 0;
}
