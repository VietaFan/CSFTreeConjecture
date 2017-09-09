#ifndef _PRINTERS_H
#define _PRINTERS_H
#include <ostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

template<typename T>
std::ostream& operator<<(std::ostream &out, std::vector<T> &vec);

template<typename T>
std::ostream& operator<<(std::ostream &out, std::vector<T> &&vec);

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream &out, std::unordered_map<T1, T2> &dict);

template<typename T>
std::ostream& operator<<(std::ostream &out, std::unordered_set<T> &S);

#endif
