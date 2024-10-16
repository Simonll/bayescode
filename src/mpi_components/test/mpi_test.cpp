#include <sstream>
#include "mpi_components/PartitionedBufferManager.hpp"
#include "mpi_components/broadcast.hpp"
#include "mpi_components/gather.hpp"
#include "mpi_components/reduce.hpp"
#include "mpi_components/scatter.hpp"
#include "operations/proxies.hpp"

using MPI::p;
using namespace std;

// clang-format off
struct MyStruct {
    int a;
    double b;
    double c;  // unused
    template <class T> void serialization_interface(T& x) { x.add(a, b); }
};
template<> struct has_custom_serialization<MyStruct> : public true_type {};
// clang-format on

struct DummyModel {
    double a{-1}, b{-1}, c{-1};
    vector<double> v, v2;
    Partition partition{{"d", "e", "f"}, static_cast<size_t>(p->size) - 1, 1};
    MyStruct j{-1, -1, -1};
    MyStruct g{-1, -1, -1};
    double h{-1};
    MyStruct i{-1, -1, -1};

    void print() {
        stringstream ss;
        for (auto e : v) { ss << e << " "; }
        stringstream ss2;
        for (auto e : v2) { ss2 << e << " "; }
        p->message(
            "Model state is:\n\ta  =  {:.2f}\n\tb  =  {:.2f}\n\tc  =  {:.2f}\n\tv  = {{{}}}\n\tv2 "
            "= "
            "{{{}}}\n\tg  = ({}, {:.2f}, {:.2f})\n\th  =  {:.2f}\n\ti  = ({}, {:.2f})\n\tj  = ({}, "
            "{:.2f})",
            a, b, c, ss.str().c_str(), ss2.str().c_str(), g.a, g.b, g.c, h, i.a, i.b, j.a, j.b);
    }
};

// clang-format off
struct StructTheGreat {
    double a;
    double b;
    template <class T> void serialization_interface(T& x) { x.add(a, b); }
};
template<> struct has_custom_serialization<StructTheGreat> : public true_type {};
// clang-format on

void compute(int, char**) {
    if (!p->rank) {
        SendBuffer buf;
        int i[4] = {2, 3, 4, 5};
        double j[4] = {2, 3, 4, 5.2};
        buf.pack(i, 4);
        buf.pack(j, 4);

        ReceiveBuffer rcvbuf(buf.size());
        memcpy(rcvbuf.data(), buf.data(), buf.size());
        p->message("Unpacked ints: {}, {}, {}, {}", rcvbuf.unpack<int>(), rcvbuf.unpack<int>(),
            rcvbuf.unpack<int>(), rcvbuf.unpack<int>());
        p->message("Unpacked doubles: {:.2f}, {:.2f}, {:.2f}, {:.2f}", rcvbuf.unpack<double>(),
            rcvbuf.unpack<double>(), rcvbuf.unpack<double>(), rcvbuf.unpack<double>());

        ReceiveBuffer rcvbuf2(buf.size());
        memcpy(rcvbuf2.data(), buf.data(), buf.size());
        auto v1 = rcvbuf2.unpack_vector<int>(4);
        p->message("Unpacked ints: {}, {}, {}, {}", v1.at(0), v1.at(1), v1.at(2), v1.at(3));
        auto v2 = rcvbuf2.unpack_vector<double>(4);
        p->message("Unpacked doubles: {:.2f}, {:.2f}, {:.2f}, {:.2f}", v2.at(0), v2.at(1), v2.at(2),
            v2.at(3));
    }
    if (!p->rank) {
        int i{1}, j{2}, k{3};
        double l{2.1}, m{3.2}, n{5.6};
        vector<int> v{5, 8, 9};
        vector<double> v2{2.35, 5.68};

        BufferManager b;
        b.add(i, j, k, l, m, n, v, v2);

        void* buf = b.send_buffer();
        size_t buf_size = b.buffer_size();
        p->message("Size of buffer is {}", buf_size);

        void* rcvbuf = b.receive_buffer();
        memcpy(rcvbuf, buf, buf_size);

        b.receive();
        p->message("{}, {}, {}, {:.2f}, {:.2f}, {:.2f}, {{{}, {}, {}}}, {{{:.2f}, {:.2f}}}", i, j,
            k, l, m, n, v.at(0), v.at(1), v.at(2), v2.at(0), v2.at(1));
    }
    if (!p->rank) {
        IndexSet is{"a", "b", "c"};
        Partition part(is, 2);
        PartitionedBufferManager bm(part);

        vector<vector<double>> v = {{1.1, 1.11}, {2.2, 2.22}, {3.3, 3.33}};
        bm.add(v);

        vector<StructTheGreat> v2 = {{0.01, 0.1}, {0.02, 0.2}, {0.03, 0.3}};
        bm.add(v2);

        auto void_buf = bm.send_buffer();
        auto buf = static_cast<double*>(void_buf);
        p->message(
            "Buf = {:.2f}, {:.2f}, {:.2f}, {:.2f} | {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}, "
            "{:.2f}, {:.2f}, {:.2f}",
            buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7], buf[8], buf[9], buf[10],
            buf[11]);
    }

    DummyModel m;
    if (!p->rank) {  // master
        m.a = 2.2;
        m.b = 3.3;
        m.c = 4.4;
        m.i = {2, 3.2, -1};
        m.j = {7, 2.21, -1};
        m.v = {-1, -1, -1};
        m.v2 = {13, 14, 15};
    } else {  // slave
        m.v = vector<double>(m.partition.my_partition_size(), p->rank + 1.1);
        m.g = {p->rank + 1, p->rank - 0.72, p->rank + 0.4};
        m.h = p->rank + 2.3;
        m.v2 = vector<double>(m.partition.my_partition_size(), -1);
    }

    Group master_operations{reduce(m.g, m.h), gather(m.partition, m.v)};

    Group slave_operations{broadcast(m.a, m.c), broadcast(m.i, m.j), scatter(m.partition, m.v2)};

    p->rank ? slave_operations.acquire() : slave_operations.release();
    p->rank ? master_operations.release() : master_operations.acquire();

    m.print();
}

int main(int argc, char** argv) { mpi_run(argc, argv, compute); }