#ifndef AMGCL_PRECONDITIONER_DUMMY2_HPP
#define AMGCL_PRECONDITIONER_DUMMY2_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   amgcl/preconditioner/dummy2.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Dummy2 preconditioner (identity matrix).
 */

#include <memory>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
    namespace preconditioner {

        template <class Backend>
        class dummy2 {
        public:
            typedef Backend backend_type;

            typedef typename Backend::matrix  matrix;
            typedef typename Backend::vector  vector;
            typedef typename Backend::value_type value_type;
            typedef typename backend::builtin<value_type>::matrix build_matrix;

            typedef amgcl::detail::empty_params params;
            typedef typename Backend::params backend_params;

            template <class Matrix>
            dummy2(
                const Matrix& M = nullptr,
                const params & = params(),
                const backend_params& bprm = backend_params()
            )
                : A(nullptr)
            {
            }

            dummy2(
                std::shared_ptr<build_matrix> M = nullptr,
                const params & = params(),
                const backend_params& bprm = backend_params()
            )
                : A(nullptr)
            {
            }

            template <class Vec1, class Vec2>
            void apply(const Vec1& rhs, Vec2&& x) const {
                backend::copy(rhs, x);
            }

            std::shared_ptr<matrix> system_matrix_ptr() const {
                return A;
            }

            const matrix& system_matrix() const {
                return *A;
            }

            size_t bytes() const {
                return 0;
            }
        private:
            std::shared_ptr<matrix>   A;

            friend std::ostream& operator<<(std::ostream& os, const dummy2& p) {
                os << "identity matrix as preconditioner" << std::endl;
                os << "  unknowns: " << backend::rows(p.system_matrix()) << std::endl;
                os << "  nonzeros: " << backend::nonzeros(p.system_matrix()) << std::endl;

                return os;
            }
        };

    } // namespace preconditioner
} // namespace amgcl
#endif
