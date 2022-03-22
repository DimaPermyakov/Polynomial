// Copyright 2022 Dmitriy <dimapermyakov55@gmail.com>
#pragma once
#include <iostream>
#include <vector>
#include <cmath>

template <class T>
class Polynomial {
public:
    Polynomial();
    Polynomial(const Polynomial<T>& polynom);
    explicit Polynomial(const std::vector<T>& vec);
    ~Polynomial();
    explicit Polynomial(const size_t& size);
    size_t Size() const;
    Polynomial<T>& operator+=(const T& polinom);
    Polynomial<T>& operator+=(const Polynomial<T>& polinom);
    Polynomial<T>& operator-=(const Polynomial<T>& polynom);
    Polynomial<T>& operator-=(const T& polynom);
    Polynomial<T>& operator*=(const Polynomial<T>& polynom);
    Polynomial<T>& operator/=(const Polynomial<T>& polynom);
    Polynomial<T>& operator/=(const T& polynom);
    Polynomial<T>& operator%=(const Polynomial<T>& polynom);
    Polynomial<T>& operator*=(const T& polynom);
    Polynomial<T>& operator=(const Polynomial<T>& polinom);
    T operator[](size_t index) const;
    T& operator[](size_t index);
    void Degree(const size_t& degree);
    void MultiplyConst(const T& factor);
    T PutUnknown(const T& x);
    void SendFactor(const size_t& degree, const T& factor);
    T GetFactor(const size_t& degree);
    void Print() const;

private:
    T* data__;
    size_t size__;
};
template <class T>
Polynomial<T>& Polynomial<T>::operator=(const Polynomial<T>& polinom) {
    if (this != &polinom) {
        delete[] data__;
        data__ = new T[polinom.size__];
        size__ = polinom.size__;
        for (size_t i = 0; i < size__; ++i) data__[i] = polinom.data__[i];
    }
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& polynom) {
    *this = *this - polynom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator-=(const T& polynom) {
    *this = *this - polynom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& polynom) {
    *this = *this * polynom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator*=(const T& polynom) {
    *this = *this * polynom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& polinom) {
    *this = *this + polinom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator+=(const T& polinom) {
    *this = *this + polinom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator/=(const Polynomial<T>& polynom) {
    *this = *this / polynom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator/=(const T& polynom) {
    *this = *this / polynom;
    return *this;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator%=(const Polynomial<T>& polynom) {
    *this = *this % polynom;
    return *this;
}
template <class T>
Polynomial<T> operator/(const Polynomial<T>& p1, const Polynomial<T>& p2) {
    size_t NewSize = p1.Size() - p2.Size() + 1;
    size_t n = p1.Size();
    size_t m = p2.Size();
    Polynomial<T> help(n);
    Polynomial<T> time(p1);
    Polynomial<T> question(NewSize);
    for (size_t i = 0; i < NewSize; ++i) {
        double coef = time[i] / p2[0];
        question[i] = coef;
        size_t k1 = 0;
        for (size_t j = 0; j < m; ++j) {
            help[k1] = coef * p2[j];
            ++k1;
        }
        // time -= help;
        // Subtraction without the displacement
        size_t max;
        size_t min;
        if (time.Size() > help.Size()) {
            max = time.Size();
            min = help.Size();
        }
        else {
            max = help.Size();
            min = time.Size();
        }
        Polynomial<T> Polinom(max);
        size_t index = 0;
        size_t index_a;
        size_t index_b;
        if (time.Size() > help.Size()) {
            for (; index < max - min; ++index) Polinom[index] = time[index];
            index_a = index;
            index_b = 0;
        }
        else {
            for (; index < max - min; ++index) Polinom[index] = help[index];
            index_a = 0;
            index_b = index;
        }
        for (size_t j = 0; j < min; ++j) {
            Polinom[index] = time[index_a] - help[index_b];
            ++index_a;
            ++index_b;
            ++index;
        }
        time = Polinom;
    }
    return question;
}
template <class T>
Polynomial<T> operator%(const Polynomial<T>& p1, const Polynomial<T>& p2) {
    Polynomial<T> remainder = p1 - p2 * (p1 / p2);
    size_t size = remainder.Size();
    if (remainder[0] == 0 && size != 1) {
        while (remainder[0] == 0) {
            for (size_t i = 0; i < size; ++i) remainder[i] = remainder[i + 1];
            --size;
        }
    }
    Polynomial<T> remainderNew(size);
    for (size_t i = 0; i < size-1; ++i) remainderNew[i] = remainder[i];
    return remainderNew;
}
template <class T>
Polynomial<T> operator/(const Polynomial<T>& p1, const T& p2) {
    Polynomial<T> quotient(p1.Size());
    for (size_t i = 0; i < p1.Size(); ++i) quotient[i] = p1[i] / p2;
    return quotient;
}
template <class T>
Polynomial<T> operator+(const Polynomial<T>& a, const Polynomial<T>& b) {
    size_t max;
    size_t min;
    if (a.Size() > b.Size()) {
        max = a.Size();
        min = b.Size();
    }
    else {
        max = b.Size();
        min = a.Size();
    }
    Polynomial<T> SumPolinom(max);
    size_t i = 0;
    size_t index_a;
    size_t index_b;
    if (a.Size() > b.Size()) {
        for (; i < max - min; ++i) SumPolinom[i] = a[i];
        index_a = i;
        index_b = 0;
    }
    else {
        for (; i < max - min; ++i) SumPolinom[i] = b[i];
        index_a = 0;
        index_b = i;
    }
    for (size_t j = 0; j < min; ++j) {
        SumPolinom[i] = a[index_a] + b[index_b];
        ++index_a;
        ++index_b;
        ++i;
    }
    return SumPolinom;
}
template <class T>
Polynomial<T> operator+(const Polynomial<T>& a, const T& b) {
    Polynomial<T> SumPolinom(a);
    SumPolinom[SumPolinom.Size() - 1] += b;
    return SumPolinom;
}
template <class T>
Polynomial<T> operator-(const Polynomial<T>& a, const Polynomial<T>& b) {
    size_t max;
    size_t min;
    if (a.Size() > b.Size()) {
        max = a.Size();
        min = b.Size();
    }
    else {
        max = b.Size();
        min = a.Size();
    }
    Polynomial<T> Polinom(max);
    size_t i = 0;
    size_t index_a;
    size_t index_b;
    if (a.Size() > b.Size()) {
        for (; i < max - min; ++i) Polinom[i] = a[i];
        index_a = i;
        index_b = 0;
    }
    else {
        for (; i < max - min; ++i) Polinom[i] = b[i];
        index_a = 0;
        index_b = i;
    }
    for (size_t j = 0; j < min; ++j) {
        Polinom[i] = a[index_a] - b[index_b];
        ++index_a;
        ++index_b;
        ++i;
    }
    size_t size = Polinom.Size();
    if (Polinom[0] == 0 && size != 1) {
        while (Polinom[0] == 0) {
            for (size_t j = 0; j < size-1; ++j) Polinom[j] = Polinom[j + 1];
            --size;
        }
    }
    Polynomial<T> PolinomNew(size);
    for (size_t j = 0; j < size; ++j) PolinomNew[j] = Polinom[j];
    return PolinomNew;
}
template <class T>
Polynomial<T> operator-(const Polynomial<T>& a, const T& b) {
    Polynomial<T> Polinom(a);
    Polinom[Polinom.Size() - 1] -= b;
    return Polinom;
}
template <class T>
Polynomial<T> operator*(const Polynomial<T>& polynom1,
    const Polynomial<T>& polynom2) {
    Polynomial<T> MultiplyPol(polynom1.Size() + polynom2.Size() - 1);
    for (size_t i = 0; i < polynom1.Size(); ++i) {
        for (size_t j = 0; j < polynom2.Size(); ++j)
            MultiplyPol[i + j] += polynom1[i] * polynom2[j];
    }
    return MultiplyPol;
}
template <class T>
Polynomial<T> operator*(const Polynomial<T>& polynom, const T& value) {
    Polynomial<T> MultiplyPol(polynom);
    MultiplyPol.MultiplyConst(value);
    return MultiplyPol;
}
template <class T>
T& Polynomial<T>::operator[](size_t index) {
    return data__[index];
}
template <class T>
T Polynomial<T>::operator[](size_t index) const {
    return data__[index];
}
template <class T>
Polynomial<T>::Polynomial(const size_t& size) : size__(size) {
    data__ = new T[size__];
    for (size_t i = 0; i < size__; ++i) data__[i] = 0;
}
template <class T>
Polynomial<T>::Polynomial(const std::vector<T>& vec) : size__(vec.size()) {
    data__ = new T[size__];
    for (size_t i = 0; i < size__; i++) data__[i] = vec[i];
}
template <class T>
Polynomial<T>::Polynomial(const Polynomial<T>& polynom)
    : size__(polynom.size__) {
    data__ = new T[size__];
    for (size_t i = 0; i < polynom.size__; ++i) data__[i] = polynom.data__[i];
}
template <class T>
Polynomial<T>::Polynomial() : size__(1) {
    data__ = new T[size__];
    data__[0] = 0;
}
template <class T>
Polynomial<T>::~Polynomial() {
    delete[] data__;
}
template <class T>
size_t Polynomial<T>::Size() const {
    return size__;
}
template <class T>
T Polynomial<T>::GetFactor(const size_t& degree) {
    return data__[size__ - 1 - degree];
}
template <class T>
void Polynomial<T>::SendFactor(const size_t& degree, const T& factor) {
    data__[size__ - 1 - degree] *= factor;
}
template <class T>
T Polynomial<T>::PutUnknown(const T& x) {
    size_t degree = size__ - 1;
    T question = 0;
    for (size_t i = 0; i < size__; ++i) {
        question += data__[i] * std::pow(x, degree);
        --degree;
    }
    return question;
}
template <class T>
void Polynomial<T>::MultiplyConst(const T& factor) {
    for (size_t i = 0; i < size__; ++i) data__[i] *= factor;
}
template <class T>
void Polynomial<T>::Degree(const size_t& degree) {
    Polynomial<T> help(size__);
    // copy the data
    for (size_t i = 0; i < size__; ++i) help.data__[i] = this->data__[i];
    size_t temp_degre = degree;
    while (temp_degre != 1) {
        *this *= help;
        --temp_degre;
    }
}
template <class T>
void Polynomial<T>::Print() const {
	size_t pow = size__ - 1;
	if (pow == 0 && data__[0] == 0) {
		std::cout << "0";
		return;
	}
	for (size_t i = 0; i < size__; ++i) {
		if (data__[i] == 0) {
			--pow;
			continue;
		}
		
		if (pow == 0)
			std::cout << std::abs(data__[i]);
		else if (pow == 1) {
			if (data__[i] == 1)
				std::cout << "x";
			else
				std::cout << std::abs(data__[i]) << "x";
		}

		else {
			if (data__[i] == 1)
				std::cout << "x^" << pow;
			else
				std::cout << std::abs(data__[i]) << "x^" << pow;
		}
		--pow;

		if (i != size__ - 1) {
			if (data__[i + 1] == 0 && i + 1 == size__ - 1)
				break;
			if (data__[i+1] < 0)
				std::cout << " - ";
			else if (data__[i + 1] == 0) {
				if (data__[i + 2] < 0)
					std::cout << " - ";
				else if (data__[i + 2] == 0)
					continue;
				else
					std::cout << " + ";

			}
				
			else
				std::cout << " + ";
		}
			
	}
}
