#pragma once

template<typename Type>
class Polynomial
{
public:
	Polynomial();
	Polynomial(const Polynomial& other);
	Polynomial& operator=(const Polynomial& rhs);
	~Polynomial();

	void SetMaxPow(int power);
	const int GetMaxPow() const;
	void AddCoef(Type coef, int power);
	const Type* GetCoefficients() const;

public:
	Type& operator[](int idx);
	const Type& operator[](int idx) const;
	Polynomial& operator+=(const Polynomial& p);
	Polynomial operator+(const Polynomial& p) const;
	Polynomial& operator-=(const Polynomial& p);
	Polynomial operator-(const Polynomial& p) const;
	Polynomial& operator*=(const Polynomial& p);
	Polynomial operator*(const Polynomial& p) const;
	Polynomial& operator*=(Type x);
	Polynomial operator*(Type x) const;
	Polynomial& operator/=(const Polynomial& p);
	Polynomial operator/(const Polynomial& p) const;
	Polynomial& operator%=(const Polynomial& p);
	Polynomial operator%(const Polynomial& p) const;
	Polynomial& operator/=(Type x);
	Polynomial operator/(Type x) const;
	Type operator()(Type x);
	Type operator()(Type down, Type up);
	Polynomial& operator++();
	Polynomial operator++(int);
	Polynomial& operator--();
	Polynomial operator--(int);
	explicit operator bool();
	bool operator !();
	explicit operator int();

	class Iterator
	{
	public:
		Iterator(Type* n) : iter(n) {}

		Type& operator*()
		{
			return *iter;
		}

		Type* operator->()
		{
			return iter;
		}

		Iterator& operator++()
		{
			++iter;
			return *this;
		}

		Iterator operator++(int)
		{
			Iterator nIter(*this);
			++(*this);
			return nIter;
		}

		Iterator& operator--()
		{
			--iter;
			return *this;
		}

		Iterator operator--(int)
		{
			Iterator nIter(*this);
			--(*this);
			return nIter;
		}

		friend bool operator==(const Iterator& lhs, const Iterator& rhs)
		{
			return lhs.iter == rhs.iter;
		}

		friend bool operator!=(const Iterator& lhs, const Iterator& rhs)
		{
			return !(lhs == rhs);
		}

		friend bool operator<(const Iterator& lhs, const Iterator& rhs)
		{
			return lhs.iter < rhs.iter;
		}

		friend bool operator>(const Iterator& lhs, const Iterator& rhs)
		{
			return !(lhs.iter < rhs.iter);
		}

		friend bool operator<=(const Iterator& lhs, const Iterator& rhs)
		{
			return !(lhs.iter > rhs.iter);
		}

		friend bool operator>=(const Iterator& lhs, const Iterator& rhs)
		{
			return !(lhs.iter < rhs.iter);
		}

	private:
		friend class Polynomial;
		Type* iter;
	};

	Iterator begin()
	{
		return Iterator(coefficients);
	}

	Iterator end()
	{
		return Iterator(coefficients + maxPow + 1);
	}

private:
	void copyFrom(const Polynomial& other);
	void free();
	void resizeToMaxPow(int newMaxPow);
	void reduceToMaxPow(int newMaxPow, int fromIdx);

private:
	Type* coefficients; // the first element is the leading coefficient
						// and the last element of the array is the coefficient for x^0 
	int maxPow;
};

template<typename Type>
Polynomial<Type>::Polynomial(): maxPow(0), coefficients(nullptr)
{
	coefficients = new(std::nothrow) Type[1];
	if (!coefficients)
	{
		std::cerr << "Not enough memory!\n";
		exit(1);
	}
	coefficients[0] = 0;
}

template<typename Type>
Polynomial<Type>::Polynomial(const Polynomial<Type>& other) : coefficients(nullptr)
{
	copyFrom(other);
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator=(const Polynomial<Type> & rhs)
{
	if (this != &rhs)
	{
		free();
		copyFrom(rhs);
	}
	return *this;
}

template<typename Type>
Polynomial<Type>::~Polynomial()
{
	free();
}

template<typename Type>
void Polynomial<Type>::SetMaxPow(int power)
{
	if (power < 0 || power < maxPow)
		return;

	maxPow = power;
}

template<typename Type>
const int Polynomial<Type>::GetMaxPow() const
{
	return maxPow;
}

template<typename Type>
void Polynomial<Type>::AddCoef(Type coef, int power)
{
	if (maxPow < power) // if the array doesn't have allocated memory for x^power's coef
	{
		resizeToMaxPow(power); // allocating memory for the coefficient
	}

	int idx = -1; // -1 because we have x^0 so we have to conform the indexes in the array 
	for (int i = 0; i <= power; ++i)
	{
		idx++;
	}

	coefficients[maxPow - idx] = coef; // maxPow - idx because of the indexation
}

template<typename Type>
const Type* Polynomial<Type>::GetCoefficients() const
{
	return coefficients;
}

template<typename Type>
Type& Polynomial<Type>::operator[](int idx)
{
	if (idx > maxPow || idx < 0)
		throw std::out_of_range("No such pow!\n");

	return coefficients[maxPow - idx];
}

template<typename Type>
const Type& Polynomial<Type>::operator[](int idx) const
{
	if (idx > maxPow || idx < 0)
		throw std::out_of_range("No such pow!\n");

	return coefficients[maxPow - idx];
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator+=(const Polynomial& p)
{
	if (maxPow == p.GetMaxPow()) // if the arrays are with the same length
	{
		for (int i = 0; i <= maxPow; ++i)
		{
			coefficients[i] += p.GetCoefficients()[i];
		}
		return *this;
	}
	else if (maxPow > p.GetMaxPow()) // if lhs is longer
	{
		int j = p.GetMaxPow() - 1;  // where the index of the leading coefficient should be
		for (int i = 0; i <= p.GetMaxPow(); ++i)
		{
			coefficients[j] += p.GetCoefficients()[i];
			j++;
		}
		return *this;
	}
	else if (maxPow < p.GetMaxPow()) //if rhs is longer
	{
		int newPow = p.GetMaxPow(); 
		Type* newCoefs = new(std::nothrow) Type[newPow + 1]; // allocating memory for coefs
		if (!newCoefs)
		{
			std::cerr << "Could not allocate memory!\n";
			exit(1);
		}
		for (int k = 0; k <= newPow; ++k) // setting coefficients to 0
		{
			newCoefs[k] = 0;
		}

		for (int i = 0; i <= p.GetMaxPow(); ++i)
		{
			newCoefs[i] = p.GetCoefficients()[i]; // adding "other's" array
		}

		int thisMaxPow = maxPow - 1; // from where the coefficients should be placed
		for (int i = 0; i <= maxPow; ++i)
		{
			newCoefs[thisMaxPow] += coefficients[i]; // from thisMaxPow to the end
			thisMaxPow++;
		}

		std::cout << std::endl;
		free();					 // deallocating
		SetMaxPow(newPow);		 // setting new maxPow
		coefficients = newCoefs; // readirecting pointer 
	}
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator+(const Polynomial& p) const
{
	Polynomial res(*this);
	res += p;
	return res;
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator-=(const Polynomial& p)
{
	if (p.GetMaxPow() < maxPow) // if rhs is longer => no problem
	{
		int pMaxPow = p.GetMaxPow() - 1; // setting index 
		for (int i = 0; i <= p.GetMaxPow(); ++i)
		{
			this->coefficients[pMaxPow] -= p.coefficients[i];
			pMaxPow++;
		}
	}
	else if (maxPow < p.GetMaxPow()) // if lhs is shorter
	{
		int tmp = p.GetMaxPow(); 
		int newPow = p.GetMaxPow(); // the new power of the polynomial
		Type* res = new(std::nothrow) Type[newPow + 1];
		if (!res)
		{
			std::cerr << "Could not allocate memory!\n";
			exit(1);
		}
		for (int i = 0; i <= newPow; ++i) // setting coeffs to 0
		{
			res[i] = 0;
		}

		int j = 0;
		for (int i = 0; i <= newPow; ++i) 
		{
			if (maxPow < tmp) // if it's rhs
			{
				res[i] -= p.GetCoefficients()[i];
				tmp--;
			}
			else // if it's lhs
			{
				res[i] = coefficients[j] - p.GetCoefficients()[i];
				j++;
			}
		}
		free();
		maxPow = newPow;
		coefficients = res;
	}
	else // if lhs and rhs are with the same length
	{
		for (int i = 0; i <= p.GetMaxPow(); ++i)
		{
			coefficients[i] -= p.GetCoefficients()[i]; // substracting
		}

		bool flag = false; // flag to see if we have to reduce the power of the polynomial
		if (coefficients[0] == 0) // if the leading coefficient is 0
			flag = true;

		int idx = 0;
		for (int j = 0; j <= maxPow; ++j)
		{
			if (coefficients[idx] == 0 && flag == true)
			{
				reduceToMaxPow(maxPow - 1, 1); // reducing to power-1 from the leading coef
				flag = true;
				continue;
			}
			flag = false;
		}

		if (coefficients[0] == 0)
		{
			reduceToMaxPow(maxPow - 1, 1);
		}
	}
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator-(const Polynomial& p) const
{
	Polynomial res(*this);
	res -= p;
	return res;
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator*=(const Polynomial& p)
{
	int newMaxPow = this->maxPow + p.GetMaxPow(); // the new pow is sum of the pows
	Type* res = new(std::nothrow) Type[newMaxPow + 1];
	if (!res)
	{
		std::cerr << "Not enough memory to multiply the Polynomials!\n";
		exit(1);
	}
	for (int c = 0; c <= newMaxPow; ++c) // setting coeffs to 0
	{
		res[c] = 0;
	}

	for (int i = 0; i <= maxPow; ++i)
	{
		for (int j = 0; j <= p.maxPow; ++j)
		{
			res[i + j] += this->coefficients[i] * p.coefficients[j]; // multiplicating
		}
	}

	free(); // deallocating
	maxPow = newMaxPow; // setting the new power
	coefficients = res; // redirecting pointer
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator*(const Polynomial& p) const
{
	Polynomial res(*this);
	res *= p;
	return res;
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator*=(Type x)
{
	for (int i = 0; i <= maxPow; ++i)
	{
		coefficients[i] *= x; // multiplicating coeffs with "x"
	}
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator*(Type x) const
{
	Polynomial res(*this);
	res *= x;
	return res;
}

template <typename Type>
Polynomial<Type>& Polynomial<Type>::operator/=(const Polynomial& p)
{
	if (this->maxPow < p.GetMaxPow())
	{
		std::cout << "Pow of divisor > pow of dividend!\n";
		return *this;
	}

	int newPow = maxPow - p.GetMaxPow(); // the new pow is difference of the pows
	Type coef; // the coefficient that we add in the quotient
	Polynomial quotient;
	int j = 0,
		i = 0;

	while(this->GetMaxPow() >= p.GetMaxPow()) // residue's pow is >= from divisor's pow
	{
		Polynomial temp; // here we add the coef^newPow 
		coef = coefficients[j] / p.GetCoefficients()[j]; // calculating coefficient
		quotient.AddCoef(coef, newPow - i); // adding the coefficient in the right place
		temp.AddCoef(coef, newPow - i); // adding the coef
		*this -= (temp*p); // multiplicating the coef^n by "p" and substracting it from
		i++;			   // this so that we can calculate the residue
	}
	*this = quotient;
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator/(const Polynomial& p) const
{
	Polynomial res(*this);
	res /= p;
	return res;
}

template<typename Type>
Polynomial<Type> & Polynomial<Type>::operator%=(const Polynomial& p)
{
	int newPow = maxPow - p.GetMaxPow(); // calculating power
	Type coef;
	int j = 0,
		i = 0;

	while (this->GetMaxPow() >= p.GetMaxPow())
	{
		Polynomial temp;
		coef = coefficients[j] / p.GetCoefficients()[j]; // calculating coef
		temp.AddCoef(coef, newPow - i); // adding the coef in front of the right power
		*this -= (temp*p); // calculating the residue
		i++;
	}
	
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator%(const Polynomial& p) const
{
	Polynomial res(*this);
	res %= p;
	return res;
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator/=(Type x)
{
	for (int i = 0; i <= maxPow; ++i)
	{
		coefficients[i] /= x; // dividing the coefficients by x
	}
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator/(Type x) const
{
	Polynomial res(*this);
	res /= x;
	return res;
}

template<typename Type>
Type Polynomial<Type>::operator()(Type x)
{
	Type value = 0; // the value in the point x
	int power; // helper for x^n = x*x*...*x (n times)
	Type tempX; // the value of x powered
	for (int i = 0; i <= maxPow; ++i)
	{
		tempX = 1; // 1 because we multiplicate and 1 is the neutral element
		power = maxPow - i; // setting power
		while (power > 0)
		{
			tempX *= x; // powering x
			power--;
		}
		value += coefficients[i] * tempX; // calculating the value in the point
	}
	return value;
}

template<typename Type>
Type Polynomial<Type>::operator()(Type down, Type up)
{
	Type downVal; // from
	Type upVal; // to
	Polynomial newPoly(*this); // the polynomial after we calculate the itegral 
	++newPoly; // calculating the integral
	upVal = newPoly(up); // calculating integral in "to" point
	downVal = newPoly(down); // calculating integral in "from" point
	return upVal - downVal;
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator++()
{
	Type coef; // the new coef
	int power; // the new power
	Polynomial newPoly; // the polynomial after the caluclation of the integral
	for (int i = 0; i <= maxPow; ++i)
	{
		if (coefficients[i] == 0)
			continue;

		power = maxPow - i + 1; // the new power is power + 1
		coef = coefficients[i] / power; // the new coef is coef / (old)power+1
		newPoly.AddCoef(coef, power); // adding the coefficient at the right power
	}
	*this = newPoly;
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator++(int)
{
	Polynomial res(*this);
	++(*this);
	return res;
}

template<typename Type>
Polynomial<Type>& Polynomial<Type>::operator--()
{
	Polynomial newPolyn; // the polynomial after we calculate derivative
	Type coef; // the new coefficients ("old" coeffs multiplicated by the "old" power)
	int power; // the new power (old power - 1)
	for (int i = 0; i <= maxPow; ++i)
	{
		power = maxPow - i - 1; // power - 1
		if (power < 0)
			break;

		coef = coefficients[i] * (maxPow - i); // calculating coefficients
		newPolyn.AddCoef(coef, power); // addignt he coefficients at the right power
	}
	*this = newPolyn;
	return *this;
}

template<typename Type>
Polynomial<Type> Polynomial<Type>::operator--(int)
{
	Polynomial res(*this);
	--(*this);
	return res;
}

template<typename Type>
Polynomial<Type>::operator bool()
{
	return maxPow == 0; // it's the zero polynomial
}

template<typename Type>
bool Polynomial<Type>::operator!()
{
	return maxPow != 0; // it's not the zero polynomial
}

template<typename Type>
Polynomial<Type>::operator int()
{
	return maxPow; // the power of the polynomial
}

template <typename Type>
void Polynomial<Type>::copyFrom(const Polynomial<Type>& other)
{
	coefficients = new(std::nothrow) Type[other.maxPow + 1]; // allocating memory
	if (!coefficients)
	{
		std::cerr << "Could not allocate memory!\n";
		exit(1);
	}

	maxPow = other.maxPow; // setting maxPow

	for (int i = 0; i <= maxPow; ++i)
	{
		coefficients[i] = other.coefficients[i]; // setting coefficients
	}
}

template<typename Type>
void Polynomial<Type>::free()
{
	delete[] coefficients;
	maxPow = 0;
}

template<typename Type>
void Polynomial<Type>::resizeToMaxPow(int newMaxPow)
{
	Type* tmp = new(std::nothrow) Type[newMaxPow + 1]; // +1 because x^0
	if (!tmp)
	{
		std::cerr << "Not enough memory!\n";
		exit(1);
	}
	
	for (int i = 0; i <= newMaxPow; ++i) // setting coeffs to 0
	{
		tmp[i] = 0;
	}

	int j = newMaxPow; // the new power
	for (int i = maxPow; i >= 0; --i) // 
	{
		tmp[j] = coefficients[i]; // assigning coefficients at the right place
		j--;
	}

	maxPow = newMaxPow; // setting the new power
	delete[] coefficients; // deallocating
	coefficients = tmp; // redirecting pointer
}

template<typename Type>
void Polynomial<Type>::reduceToMaxPow(int newMaxPow, int fromIdx)
{
	Type* tmp = new(std::nothrow) Type[newMaxPow + 1]; // +1 because x^0
	if (!tmp)
	{
		std::cerr << "Not enough memory!\n";
		exit(1);
	}

	for (int i = 0; i <= newMaxPow; ++i) // setting coefficients to 0
	{
		tmp[i] = 0;
	}

	for (int i = 0; i <= newMaxPow; ++i)
	{
		tmp[i] = coefficients[fromIdx]; // reducing coefficients from index
		fromIdx++;
	}

	maxPow = newMaxPow; // setting the new power
	delete[] coefficients; // deallocating
	coefficients = tmp; // redirecting pointer
}
///////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
bool operator==(const Polynomial<Type>& a, const Polynomial<Type>& b)
{
	if (a.GetMaxPow() != b.GetMaxPow()) // different powers
		return false;

	for (int i = 0; i <= a.GetMaxPow(); ++i) // comparing the coefficients
	{
		if (a[i] != b[i])
		{
			return false;
		}
	}
	return true;
}

template<typename Type>
bool operator!=(const Polynomial<Type>& a, const Polynomial<Type>& b)
{
	return !(a == b);
}

template<typename Type>
bool operator<(const Polynomial<Type>& a, const Polynomial<Type>& b)
{
	return a.GetMaxPow() < b.GetMaxPow(); // comparing the powers
}

template<typename Type>
bool operator>(const Polynomial<Type>& a, const Polynomial<Type>& b)
{
	return b.GetMaxPow() < a.GetMaxPow();
}

template<typename Type> 
bool operator<=(const Polynomial<Type>& a, const Polynomial<Type>& b)
{
	return !(a > b);
}

template<typename Type>
bool operator>=(const Polynomial<Type>& a, const Polynomial<Type>& b)
{
	return !(a < b);
}

template<typename Type>
Polynomial<Type> operator*(Type x, const Polynomial<Type>& p) // i.e. 5*polynomial
{
	Polynomial<Type> tmp(p); // copy
	for (int i = 0; i <= p.GetMaxPow(); ++i)
	{
		tmp.AddCoef(p.GetCoefficients()[i] * x, i); // adding coefficients at the pows
	}
	return tmp;
}

template<typename Type>
std::ostream& operator<<(std::ostream& out, const Polynomial<Type>& h)
{
	out << "f(x) = ";
	for (int i = 0; i <= h.GetMaxPow(); ++i) // f(x) = A0*x^n + ... + An-1*x + An
	{
		out << h.GetCoefficients()[i] << ".x^" << h.GetMaxPow() - i;
		if (i != h.GetMaxPow())
		{
			std::cout << " + ";
		}
	}
	return out;
}

template<typename Type>
std::istream& operator>>(std::istream& in, Polynomial<Type>& h)
{
	int coef;
	int pow;
	std::cout << "Enter the pow of the Polynomial: ";
	in >> pow;
	std::cout << "Enter the coefficients of the Polynomial(starting from the leading): ";
	for (int i = 0; i <= pow; ++i)
	{
		in >> coef;
		h.AddCoef(coef, pow - i);
	}
	return in;
}
