# BIE-PA2_second
### Second homework on BIE-PA2 course at FIT CTU

The task is to develop a C++ class CPolynomial representing polynomials of arbitrary degree.

The class will store individual coefficients (floating point numbers - double precision), moreover, public interface will offer overloaded operators and methods to manipulate polynomials:

default constructor
the constructor will prepare an object representing polynomial 0,
copy constructor
the constructor will be implemented if your internal representation requires it,
destructor
destructor will be present if your internal representation requires it,
overloaded operator=
the operator will copy the contents from one instance to another (if the automatically generated op= cannot be used with your implementation).
operator <<
will output a textual form of the polynomial into an output stream. The output formatting must follow these rules:
the polynomial will be displayed in a decreasing order of powers, i.e. from the highest power of x,
terms with zero coefficient are not displayed in the output,
terms with coefficient +1 and -1 are displayed without the coefficient (just a power of x),
there are not unnecessary - in the output (i.e. x - 9 or - x^2 + 4 is OK, whereas x + (-9) is not),
zero polynomial will be displayed as 0.
operator *
multiplies given polynomial either with a double, or with another polynomial,
operator *=
multiplies given polynomial either with a double, or with another polynomial,
operators == and !=
compare two polynomials for exact match (no epsilon tolerance),
operator []
is used to access (read / write) individual coefficients, the term is determined by the index (0 = absolute, 1 = x, 2 = x^2, ... ). The read variant of [] must be available even on const instances,
operator ()
evaluates the value of the polynomial for the given value x (x is a double),
degree()
the method returns the degree of the polynomial (e.g. x^3+4 has degree of 3, 5 has degree 0, specifically 0 has degree 0).
cast to bool type
returns true if the polynomial is not zero (some coefficient is not zero).
operator !
returns true if the polynomial is zero (all coefficietns are zeroes).
Optional manipulator poly_var:

the manipulator sets the name of the variable for the output. The default variable in the output is x, the manipulator may be used to set the name to an arbitrary string,
the implementation is optional. If you decide not to implement the manipulator, keep the dummy implementation from the attached example. The dummy implementation does not do anything useful, however, it is needed to pass the compilation. The compilation fails if the manipulator is not present at all,
the attached examples demonstrate how the manipulator works.
Submit a source file with your implementation of CPolynomial class (and the manipulator, if implemented). The submitted file shall not contain any #include directives nor main function. If your main function or #include remains in the file, please place them into a conditional compile block.

This task does not provide the required class interface. Instead, you are expected to develop the interface yourself. Use the description above, the attached examples, and your knowledge of overloaded operators.
