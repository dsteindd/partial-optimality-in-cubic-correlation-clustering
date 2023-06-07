#include <stdexcept>
#include <iostream>

#include "andres/partition.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed.");
}

int main() {
    typedef andres::Partition<int> Partition;

    {
        const size_t numberOfElements = 5;

        Partition partition(numberOfElements);


        test(partition.numberOfSets() == 5);
        test(partition.numberOfElements() == 5);

        partition.merge(1, 2);
        test(partition.numberOfSets() == 4);
        test(partition.numberOfElements() == 5);

        std::vector<int> r(static_cast<std::size_t>(partition.numberOfElements()));

        partition.elementLabeling(r.begin());

        std::vector<int> expected({0, 1, 1,2,3});

        int index = 0;
        for (auto l : r){
            std::cout << l << std::endl;
            test(l == expected[index]);
            index++;
        }


        partition.merge(1, 3);

        test(partition.numberOfSets() == 3);

        std::vector<int> expected2({0, 1, 1,1,2});

        std::vector<int> r2(static_cast<std::size_t>(partition.numberOfElements()));

        partition.elementLabeling(r2.begin());

        int index2 = 0;
        for (auto l : r2){
            std::cout << l << std::endl;
            test(l == expected2[index2]);
            index2++;
        }
    }
    return 0;
}