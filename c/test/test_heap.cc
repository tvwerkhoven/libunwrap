/*
 test_heap.cc -- Test heap construct using priority_queue
 Copyright (C) 2012 Visa Korkiakoski <korkiakoski@strw.leidenuniv.nl>
 & Tim van Werkhoven <werkhoven@strw.leidenuniv.nl>
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdio.h>

#include <queue>

using namespace std;

class QualityComparison {
  const double *const quality;
public:
  QualityComparison(const double *const qualitypo) : quality(qualitypo) { }
  
  bool operator() (const int& lhs, const int&rhs) const {
    return (quality[lhs] > quality[rhs]);
  }
};

class mycomparison
{
  bool reverse;
public:
  mycomparison(const bool& revparam=false)
  {reverse=revparam;}
  bool operator() (const int& lhs, const int&rhs) const
  {
    if (reverse) return (lhs>rhs);
    else return (lhs<rhs);
  }
};


int main() {
  double qual[] = {10,20,30,40,50,60,70,80,90,100};
  
  // priority_queue< int, vector<int>, QualityComparison> Adjoint;
  typedef priority_queue< int, vector<int>, QualityComparison > itohpq_t;
  itohpq_t Adjoint(qual);
  
  // using mycomparison:
  priority_queue< int, vector<int>, mycomparison > fourth(true);
  
  typedef priority_queue<int,vector<int>,mycomparison> mypq_type;
  mypq_type fifth ();
  mypq_type sixth (true);
  
  fourth.push(10);
  fourth.push(20);
  fourth.push(1);

  Adjoint.push(7);
  fprintf(stderr, "Top: %d\n", Adjoint.top());
  Adjoint.push(8);
  fprintf(stderr, "Top: %d\n", Adjoint.top());
  Adjoint.push(9);
  fprintf(stderr, "Top: %d\n", Adjoint.top());
  
  
  // Pop values
  fprintf(stderr, "pop & top: %d\n", Adjoint.top());
  Adjoint.pop();
  fprintf(stderr, "pop & top: %d\n", Adjoint.top());
  Adjoint.pop();
  fprintf(stderr, "pop & top: %d\n", Adjoint.top());
  Adjoint.pop();
  fprintf(stderr, "pop & top: %d\n", Adjoint.top());

  return 0;
}