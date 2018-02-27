//
// Created by Govinda Kamath on 2/12/18.
//

#include <boost/heap/priority_queue.hpp>
#include <vector>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>

using namespace boost::heap;


int main()
{
    fibonacci_heap<int> fh;
    auto handle = fh.push(2);
    fh.push(3);
    fh.push(1);

    std::cout << fh.top() << '\n';
    fh.update(handle, 4);

    std::cout << fh.top() << '\n';
}