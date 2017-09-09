#include "printers.h"
#include <ostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
using namespace std;

template<typename T>
ostream& operator<<(ostream &out, vector<T> &vec) {
	out << "[";
	if (vec.size() > 0) {		
		for (int i=0; i<vec.size()-1; ++i)
			out << vec[i] << ", ";
		out << vec[vec.size()-1];
	}
	out << "]";
	return out;
}

template<typename T>
ostream& operator<<(ostream &out, vector<T> &&vec) {
	out << "[";
	if (vec.size() > 0) {		
		for (int i=0; i<vec.size()-1; ++i)
			out << vec[i] << ", ";
		out << vec[vec.size()-1];
	}
	out << "]";
	return out;
}

template<typename T1, typename T2>
ostream& operator<<(ostream &out, unordered_map<T1, T2> &dict) {
	out << "{";
	int k = dict.size();
	for (auto it=dict.begin(); it!=dict.end(); ++it) {
		out << (it->first) << ": " << (it->second);
		if (k > 1) {
			out << ", ";
		}
		--k;
	}
	out << "}";
	return out;
}

template<typename T>
ostream& operator<<(ostream &out, unordered_set<T> &S) {
	out << "{";
	int k = S.size();
	for (auto it=S.begin(); it!=S.end(); ++it) {
		out << (*it);
		if (k > 1) {
			out << ", ";
		}
		--k;
	}
	out << "}";
	return out;
}

