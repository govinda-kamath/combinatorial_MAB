//
// Created by Govinda Kamath on 2/12/18.
//

#include <boost/heap/priority_queue.hpp>
#include <vector>
#include <boost/heap/binomial_heap.hpp>
#include <iostream>

using namespace boost::heap;

int main()
{
    binomial_heap<int> bh;
    auto handle = bh.push(2);
    bh.push(3);
    bh.push(1);

    bh.update(handle, 4);

    std::cout << bh.top() << '\n';
}