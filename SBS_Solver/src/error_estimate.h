#pragma once
template<class T>
class ErrorEstimate {
	T x0, x1;
public:
	ErrorEstimate();
	ErrorEstimate(const T &);
	void Init(const T &);
	void Update(const T &);
	double RelativeErr(void) const;
	static double GetRelativeMax(std::initializer_list<ErrorEstimate<T>>);
};

template<class T>
ErrorEstimate<T>::ErrorEstimate()
{

}

template<class T>
ErrorEstimate<T>::ErrorEstimate(const T& x0_)
{
	Init(x0_);
}

template<class T>
void ErrorEstimate<T>::Update(const T& x1_)
{
	x1 = x1_;
}

template<class T>
double ErrorEstimate<T>::RelativeErr(void) const
{
	if (x0 == 0.0) {
		if (x1 == 0.0)
			return (0.0);
		else
			return(1.0);
	};
	return std::abs(x0 - x1) / std::abs(x0);
}

template<class T>
double ErrorEstimate<T>::GetRelativeMax(std::initializer_list<ErrorEstimate<T>> list)
{
	double err = 0;
	for (auto it = list.begin(); it != list.end(); ++it) {
		double it_err = it->RelativeErr();
		if (err < it_err)
			err = it_err;
	}
	return (err);
}

template<class T>
inline void ErrorEstimate<T>::Init(const T & x0_)
{
	x0 = x0_;
}