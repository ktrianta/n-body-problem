#ifndef ARRAY_2D
#define ARRAY_2D

#include <vector>
#include <algorithm>

#ifdef DYNAMIC

template<class T>
class Array2D {
    public:
        using data_type = std::vector<T>;
        using value_type = typename data_type::value_type;
        using size_type = typename data_type::size_type;
        using reference = typename data_type::reference;
        using const_reference = typename data_type::const_reference;

        Array2D(size_type rows, size_type cols)
            : m_cols(cols), m_data(rows * cols)
        {}

        Array2D(size_type rows, size_type cols, const_reference val)
            : m_cols(cols), m_data(rows * cols, val)
        {}

        reference operator() (size_type const row, size_type const column)
        {
            return m_data[m_cols*row + column];
        }

        const_reference operator() (size_type const row, size_type const column) const
        {
            return m_data[m_cols*row + column];
        }

    private:
        size_type m_cols;
        data_type m_data;
};

#else

template<class T>
class Array2D {
    public:
        Array2D(size_t rows, size_t cols)
            : m_data(new T[rows][column_size])
        {}

        Array2D(size_t rows, size_t cols, const T& val)
            : m_data(new T[rows][column_size])
        {
            std::fill(&m_data[0][0], &m_data[0][0] + rows*column_size, val);
        }

        T& operator() (size_t const row, size_t const column)
        {
            return m_data[row][column];
        }

    private:
        static const size_t column_size = 2;

        T (*m_data)[column_size];
};

#endif // DYNAMIC

#endif // ARRAY_2D 
