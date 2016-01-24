/**
 * @author Mehmet Akif TAŞOVA <makiftasova@gmail.com>
 *
 * An implementation of Gauss–Legendre Algorithm with GMP Library 
 *
 */
#include <unistd.h>
#include <stdint.h>
#include <iostream>
#include <cmath>

#include <getopt.h>

#include <gmpxx.h>

/* EXIT CODES */
const uint8_t EXIT_ILLEGAL_ARGUMENT(1);

//const mpz_class NUMBER_OF_ITERATIONS(1e5);
//const mpf_class NUMBER_OF_ITERATIONS(1e-1234_mpf);
const uint64_t NUMBER_OF_ITERATIONS(15); //UINT64_MAX -1 ;

static inline mpf_class pow(const mpf_class op, const uint64_t pw){
	mpf_class res(1);
	for(uint64_t i = 0; i < pw; ++i){
		res = res * op;
	}
	return res;
}

static inline void printHelp(const char* argv0){
	std::cout << argv0 << " [-h|-v|-n <number>]" << std::endl;
	std::cout << "This program calculates Number PI with using Gauss–Legendre Algorithm" << std::endl;
	std::cout << "This program can accept the command line arguments given below" << std::endl;
	std::cout << "\t-h\t\tPrints this help message and exits." << std::endl;
	std::cout << "\t-v\t\tActivates verbose mode thus making program to print lots of things." << std::endl;
	std::cout << "\t-n <number>\tSets maximum iteration number for Gauss–Legendre Algorithm (default: " 
																		<< NUMBER_OF_ITERATIONS << ")" << std::endl;
}

int main(int argc, char ** argv){ 

	uint64_t numOfIterations = NUMBER_OF_ITERATIONS;
	bool isVerboseMode = false;

	int c;
	while ((c = getopt (argc, argv, "hvn:")) != -1)
    switch (c)
    {
		case 'h': //help
			printHelp(argv[0]);
			return EXIT_SUCCESS;
		break;

		case 'v': //verbose mode
			isVerboseMode = true;
			if(isVerboseMode)
				std::cout << "[INFO] Verbose mode activated.\n";
		break;

		case 'n':
			numOfIterations = atoi(optarg);
			if(isVerboseMode)
				std::cout << "[INFO] Will run for only  " << numOfIterations << " iterations.\n";
		break;
	    
		case '?':
			if (optopt == 'n')
				std::cout << "[ERR] Option -"<< (char)optopt << " requires an argument \n";
			else if (isprint (optopt))
				std::cout << "[ERR] Unknown option: `"<< (char)optopt << "`\n" ;
			else
				std::cout << "[ERR] Unknown option character: `"<< optopt << "`\n" ;
			return EXIT_ILLEGAL_ARGUMENT;
		default:
			abort ();
      }

	std::cout.precision(100);
	//std::cout.setf(std::ios::fixed, std::ios::floatfield);

	mpf_class a(1), b(1 / (sqrt(2))), t(1.0 / 4.0), p(1), est_pi(0), a_prev(1);

	for(uint64_t i = 0; i <numOfIterations; ++i){
	//while( abs(a - b) > NUMBER_OF_ITERATIONS ){

		a_prev = a;

		a = (a + b) / 2; // a[n+1] = (a[n] + b[n]) / 2

		b = sqrt(a_prev * b); // b[n+1] = sqrt(a[n] * b[n]) 

		t = t - (p * pow((a_prev - a) , 2) ); // t[n+1] = t[n] - p[n] * ((a[n] - a[n+1]) ** 2)

		p = 2 * p; // p[n+1] = 2 * p[n]

		est_pi =  pow((a + b), 2) / (4 * t);

		if( ((i % 10) == 0) && isVerboseMode){
			std::cout << "[PROG] Iteration: " << i << " Current Estimation: " << est_pi << std::endl;
		}

	}

	std::cout << "Estimated Value of PI: " << est_pi << std::endl;
	std::cout << "Number of Iterations: " << numOfIterations << std::endl;




	return EXIT_SUCCESS;
}