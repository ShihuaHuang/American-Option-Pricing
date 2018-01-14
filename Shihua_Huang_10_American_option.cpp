#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <vector>
using namespace std;

double up_factor, uptick_prob, downtick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;
double **mem1, **mem2;

double max(double a, double b) {
	return (b < a) ? a : b;
}

void initialize_option(int no_of_divisions)
{
	//for call option, set (no_of_divisions+1)*(2*k+1) array, k from 0 to no_of_divisions
	mem1 = new double*[no_of_divisions + 1];//for k to change, there are no_of_divisions+1 rows
	for (int a = 0; a <= no_of_divisions; a++)
		mem1[a] = new double[2 * a + 1];//for i to change, it depends on k and the relationship is 2*k+1, set columns to be this
	for (int a = 0; a <= no_of_divisions; a++)
	{
		for (int b = 0; b <= 2 * a; b++)
			mem1[a][b] = -1;//initialize every element to be -1
	}
	//for put option, set array the same way
	mem2 = new double*[no_of_divisions + 1];
	for (int a = 0; a <= no_of_divisions; a++)
		mem2[a] = new double[2 * a + 1];
	for(int a = 0; a <= no_of_divisions; a++)
	{
		for (int b = 0; b <= 2 * a; b++)
			mem2[a][b] = -1;
	}
}

double american_call_option(int k, int i, double current_stock_price) {
	if (mem1[k][k + i] != -1)
		return mem1[k][k + i];//use memoization to store value and call it in future use
	else
	{
		if (k == no_of_divisions)
			//figure out the bound or stop condition
			return mem1[k][k + i] = max(0.0, (current_stock_price - strike_price));
		else
			//plug the formula to calculate the new array for the first time
			return mem1[k][k + i] = max((current_stock_price - strike_price),
			(uptick_prob*american_call_option(k + 1, i + 1, current_stock_price*up_factor) + downtick_prob*american_call_option(k + 1, i - 1, current_stock_price / up_factor)
				+ (1 - uptick_prob - downtick_prob)*american_call_option(k + 1, i, current_stock_price)) / R);
	}

}

double american_put_option(int k, int i, double current_stock_price) {
	if (mem2[k][k + i] != -1)
		return mem2[k][k + i];
	else
	{
		if (k == no_of_divisions)
			return mem2[k][k + i] = max(0.0, (strike_price - current_stock_price));
		else
			return mem2[k][k + i] = max((strike_price - current_stock_price),
			(uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor) + downtick_prob*american_put_option(k + 1, i - 1, current_stock_price / up_factor)
				+ (1 - uptick_prob - downtick_prob)*american_put_option(k + 1, i, current_stock_price)) / R);
	}

}

int main(int argc, char* argv[])
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_time;

	start = std::chrono::system_clock::now();//start calculating the time

	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%lf", &risk_free_rate);
	sscanf_s(argv[4], "%lf", &volatility);
	sscanf_s(argv[5], "%lf", &initial_stock_price);
	sscanf_s(argv[6], "%lf", &strike_price);

	up_factor = exp(volatility*sqrt(2 * expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));
	uptick_prob = pow((sqrt(R) - 1 / sqrt(up_factor)) / (sqrt(up_factor) - 1 / sqrt(up_factor)), 2);
	downtick_prob = pow((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - 1 / sqrt(up_factor)), 2);

	cout << "(Memoizied) Recursive Trinomial American Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << 1 - uptick_prob - downtick_prob << endl;
	cout << "--------------------------------------" << endl;
	initialize_option(no_of_divisions);
	double call_price = american_call_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Call Option = " << call_price << endl;
	double put_price = american_put_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Put Option = " << put_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Let us verify the Put-Call Parity: S+P-C = Kexp(-r*T) for American Options" << endl;
	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " ?=? " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	if (abs(initial_stock_price + put_price - call_price - strike_price*exp(-risk_free_rate*expiration_time)) <= 1e-3)
		cout << "Looks like Put-Call Parity holds within three decimal places" << endl;
	else
		cout << "Looks like Put-Call Parity does NOT hold" << endl;
	cout << "--------------------------------------" << endl;

    end = std::chrono::system_clock::now();//stop calculating the time
	elapsed_time = end - start;//calculate the time used
	cout << "Elapsed time to complete: " << elapsed_time.count() << "s\n" << endl;//output the time
}