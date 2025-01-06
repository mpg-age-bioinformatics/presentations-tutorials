# ---
# jupyter:
#   jupytext:
#     cell_markers: '{{{,}}}'
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3.10.8
#     language: python
#     name: py3.10.8
# ---

# # Python workshop
# ## bioinformatics@age.mpg.de
# ## http://bioinformatics.age.mpg.de/
#
# &nbsp;
#
# ---
#

# ## Pre workshop instructions
#
# 1. Install docker: https://www.docker.com/products/docker-desktop
#
# 2. Pull the workshop image
# ```
# docker pull mpgagebioinformatics/python_workshop
# ```
#
# 3. Run the image and enter the container:
# ```
# mdkir -p ~/python_workshop
# docker run -p 8787:8787 -p 8888:8888 -v ~/python_workshop/:/home/mpiage --name python_workshop -it mpgagebioinformatics/python_workshop:latest
# ```
#
# 4. Start jupyterhub inside the container
# ```
# module load jupyterhub
# jupyter notebook --ip=0.0.0.0
# ```
#
# 5. Access jupyter.
# You will be shown a message like this:
# ```
#    Copy/paste this URL into your browser when you connect for the first time,
#     to login with a token:
#         http://b4b294d48b2b:8888/?token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478&token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478
# ```
# Replace the "b4b294d48b2b" string (ie. bettween the http:// and :8888/) with "localhost" (eg. http://localhost:8888/?token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478&token=01d1532b8624becc5c7b288cf0c600b7805e802349e81478 ) and paste it into your web browser.
#
#
# If you are not able to install docker you can still run the notebook by installing Python3 - https://www.python.org - and jupyter - https://jupyter.org. Afterwards make sure you have installed all required packages:
# ```
# Package            Version
# ------------------ -------
# autograd           1.2    
# backcall           0.1.0  
# bleach             2.1.3  
# Bottleneck         1.2.1  
# cycler             0.10.0 
# decorator          4.3.0  
# entrypoints        0.2.3  
# future             0.17.1 
# html5lib           1.0.1  
# ipykernel          4.8.2  
# ipython            6.4.0  
# ipython-genutils   0.2.0  
# ipywidgets         7.2.1  
# jedi               0.12.0 
# Jinja2             2.10   
# jsonschema         2.6.0  
# jupyter            1.0.0  
# jupyter-client     5.2.3  
# jupyter-console    5.2.0  
# jupyter-core       4.4.0  
# kiwisolver         1.0.1  
# lifelines          0.19.4 
# MarkupSafe         1.0    
# matplotlib         3.0.2  
# matplotlib-venn    0.11.5 
# mistune            0.8.3  
# nbconvert          5.3.1  
# nbformat           4.4.0  
# notebook           5.5.0  
# numpy              1.16.1 
# pandas             0.24.1 
# pandocfilters      1.4.2  
# parso              0.2.1  
# pexpect            4.6.0  
# pickleshare        0.7.4  
# pip                9.0.3  
# prompt-toolkit     1.0.15 
# ptyprocess         0.6.0  
# Pygments           2.2.0  
# pyparsing          2.3.1  
# python-dateutil    2.7.3  
# pytz               2018.9 
# pyzmq              17.0.0 
# qtconsole          4.3.1  
# scikit-learn       0.20.2 
# scipy              1.2.1  
# seaborn            0.9.0  
# Send2Trash         1.5.0  
# setuptools         39.0.1 
# simplegeneric      0.8.1  
# six                1.11.0 
# sklearn            0.0    
# terminado          0.8.1  
# testpath           0.3.1  
# tornado            5.0.2  
# traitlets          4.3.2  
# wcwidth            0.1.7  
# webencodings       0.5.1  
# widgetsnbextension 3.2.1  
# xlrd               1.2.0 
# ```
#
# Check your packages with `pip3 list` and install with `pip3 install <package_name>==<version> --user`.
#
# You can then download the notebook from https://github.com/mpg-age-bioinformatics/presentations-tutorials/blob/gh-pages/presentations/modules/python_workshop/python_workshop.ipynb and start your jupyter with `jupyter notebook --ip=0.0.0.0`.
#
# &nbsp;
#
#
# ---



# # Why Python?
# - Easy to learn
#
#
#  - It's free and well documented. 
#
#
#  - Popular (easy to get help on internet forums), big community.
#
#
#  - Fast !?
#
#
#  - Scripts are portable (can run on windows, mac os or linux).
#
#
#  - A lot of libraries for many applications ready to use.
#
# &nbsp;
#
# ---
#
# &nbsp;

# # Basic syntax
#
# &nbsp;
#
# - Python is a “high level” scripting language, witch means it is as close as possible to natural human language. 
#
#
# - It is possible to run commands interactively  or in scripts. 
#
#
# - On the terminal, type python: 
#
# ![1](1.png)
#
#
#
#

# ![2](2.png)
#
# - Print is a function, usually in python all the functions are called by name_function()
#
#
# - But before python 3, print function could be called without parenthesis ()
#
#
# - Inside the parenthesis, for some functions, there are parameters, called arguments
#
#
# - Inside a script it looks like this:
#
#

# ![3](3.png)
#
# - First line points to the python interpreter on the operational system
#
#
# - And now to run (on unix based system):
#

# ![4](4.png)
#

# - Using # on a script, you can write comments and they will be ignored for the execution of the program
#
#
# - Blank lines are ignored in a script
#
#
# - It is possible to write multiple statements on a single line, separated by ; 
#

# ![6](6.png)
#
# - **From now, in this workshop, we will only use interactive python withing a Jupyter notebook.**
#

# {{{
# First print on Jupyter

print("Hello, Python!")
# }}}

# ---
# # Variables
#
# - A variable is something which can change!! 
#
#
# - It is a way of referring to a memory location by a computer program.
#
#
# - It stores values, has a name (identifier)  and data type.
#
#
# - While the program is running, the variable can be accessed, and sometimes can be changed. 
#
#
# -  Python is not “strongly-typed” language. It means that the type of data storage on variables can be changed (other languages like C the variable should be declared with a data type).
#
# &nbsp;
#
# ---

# # Variables and identifiers
#
# -  Variables and identifiers are often a source of confusion. 
#
#
# - But identifiers are names of variables, they have name AND other features, like a value and  data type.
#
#
# - In addition, identifiers are not used only for variables, but also for functions, modules, packages, etc. 
#
#
# - A valid identifier is a non-empty sequence of characters of any length with:
#
#  a) The start character can be the underscore "_" or a capital or lower case letter.
#
#  b) The letters following the start character can be anything which is permitted as a start character plus the digits.
#
#  c) Just a warning for Windows-spoilt users: Identifiers are case-sensitive!
#
#  d) Python keywords are not allowed as identifier names! Ex: `and`, `as`, `assert`, `break`, `class`, `continue`, `def`, `del`, `elif`, `else`, `except`, `for`, `if`, `in`, `is`, `lambda`, `not`, `or`, `pass`, `return`, `try`, `with`.
#
# &nbsp;
#

# {{{
# Some example of variables

i = 30
J = 32.1

some_string= "string"
# }}}

# # Basic Data Types summary in Python
#
# ### Numerics:
#
# 1. int: integers, ex: `610`, `9580`
#
#
# 2. long: long integers of non-limited length (only python 2.x) <- Python 3 int is unlimited  
#
#
# 3. Floating-point, ex: `42.11`, `2.5415e-12`
#
#
# 4. Complex, ex: `x = 2 + 4i`
#
# ### Sequences:
#
# 1. `str`: Strings (sequence of characters), ex:  `"ABCD"`, `"Hello world!"`, `"C_x-aer"` 
#
#
# 2. `list`
#
#
# 3. `tuple`
#
#
# ### Boolean: 
# - `True` or `False`
#
# ### Mapping: 
# - `dict` (dictionary)

#
#

# # Numbers
#
# - int (integers): positive or negative numbers without decimal point.

a_number = 35
type(a_number)

#  - Note that there is no “ ”. If you include quotation marks:
#  
#  It is not a integer (int) but a string (str) type. 
#

a_number = "35" 
type(a_number)

# - float (floating point real values): real numbers with decimal points (ex: 23.567) and also can be in scientific notation (1.54e2 – it is the same of 1.54 x 100 and the same of 154) 
#

a_float_number=23.567
type(a_float_number)

another_float_number = 1.54e20
type(another_float_number)

another_float_number

another_float_number

#  - complex (complex numbers): not used much in Python. They are of the form a + bJ (a and b are floats) and J an imaginary  number

a_complexcomplex=23.567+0j

type(a_complexcomplex)

# ### Number type conversions

a_float_number

int(a_float_number)

a_float_number

# {{{
# modifications on-fly during workshop

a_int_number = int(a_float_number)
# }}}

a_int_number

complex(a_float_number)

str(a_float_number)

# # Hands-on:
#
# - Create a variable numeric type float
# - Convert it to a type inter 





# ---
# # Numerical operations
#
# - Unlike for strings, for numeric data types the operators `+`,`*`,`-`,`/` are arithmetic

x=55
y=30

x+y

x-y

x*y

x/y  # why 1 ??

float(x)/float(y)

# ####  Some extra operators:
#
# - % (Modulus) return the remainder of a subtraction
# -   ** (Exponent) return result of exponential calculation

x%y

# {{{
# explained:
a=55/30
b=int(a)
c=b * 30
d=55 - c

print(f"a: {a}; b: {b}; c: {c}; d: {d}")

# 55 - int(55/30)*30
# }}}

x**y

2**2

# # Hands-on
#
# - Create 3 variables with values 4, 12.4 and 30.
# - Calculate the sum of all and store it into a new variable
# - Multiple the result by 5 and assign it to another variable





# # Strings
#
#  - Strings are marked by quotes 
#
#
#  - Wrapped with single-quote ('):
#

'This is a string with single quotes'

#  - Wrapped with double-quote ("): 

"This is a string with double quotes"

#  -  Wrapped with three characters, using either single-quote or double-quote:
#

'''A String in triple quotes can extend over multiple lines, and can contain "double" quotes.'''

#  - A string in Python consists of a series or sequence of characters - letters, numbers, and special characters. 
#
# - Strings can be indexed. The first character of a string has the index 0.
#

str_1= "A string consists of characters"

str_1

str_1[0]

str_1[3]

len(str_1)

#  - Last character:
#

str_1[-1]

#  - This last one is possible because the index can be also counted from the right, using negative values:
#
#  "STRING"
# [-6], [-5], [-4], [-3], [-2], [-1]
#
# - In addition to the normal way, from the left:
#
#  "STRING"
# [0], [1], [2], [3], [4], [5]
#

# # Strings are “immutable”
#
# - Like in Java, Python strings cannot be changed.

# a wrong example (which would return an error):
# ```
# s = "Some things are immutable!"
# s[-1] = "."
# ```


# # Operations with strings
#
# - Concatenation: using the operator `+` it is possible to concatenate 2 or more strings: 
#
#
#
#

# #### Attention to the absence of space

"Hello" + "World"  #  <- *Attention to the absence of space*

#  - Repetition: using the operator  `*` it is possible to repeat n times one string:

"HelloWorld" * 3

#  - Indexing: as mentioned before, is possible to recover one specific position of the string by index: 

"HelloWorld"[0] 

#  - Slicing: used to recover substring of strings with slicing notation:
#

"HelloWorld"[2:4]  # Note that [2:4] will get the 3nd  and 4rd  position, not the 5th 

#  - Size

len("HelloWorld")

type(len("HelloWorld")) # This is a int number

#  - Split: It is possible to use a substring (or character) to split a string:

s = "Some things are immutable!"


print(s.split())  # Default is space character

print(s.split("a"))


s2 = "You can split some strings, and you can decide how to do it!"


print(s2.split(","))


print(s2.split(",")[0])


# # Hands-on
#
# - Print the 7th letter of the string "RafaelCuadrat"
# - Split the names  "Jorge,Rafael,Franziska,Daniel" on string `,`
# - Create a string using the 2 first and 2 last characters from the string "I love python"
# - Add "ing" to the end of the string created.





# ---
# # Escape sequences 
#
#  - If you want to use some special characters with the literal meaning inside a string, you need to use backslash ie. `\`  before. 
#
#  Ex:  `"The double quotation mark is \""` -> If you use only `"` without `\` it will close the string without show `"` inside the string
#
#  Ex2: `"I like to use \\"`   -> You need to use a backslash before a backslash to display `\` inside a string
#
# &nbsp;
#
#  - Another scape sequences using backslash include:
#
#      - New line:  `\n`ewline
#
#      - tabular space: `\t`   <- very important dealing with tabular files
#

# {{{
# a wrong example
# print("This is quotes mark " in a string")
# }}}

print("This is a quote mark \" in a string")

print("This is tabular t in a string")

print("This is tabular \t in a string")

# ---
# # Lists
#
#  -  It  is a versatile data type in python and can be written as comma-separated values (items) limited by [ ]. 
#
#
#  - Items in the list do not need to be of the same type! 
#

list1 = ['physics', 'chemistry', 1997, 2000]
list2 = [1, 2, 3, 4, 5 ]
list3 = ["a", "b", "c", "d"]


#  - Similar to string, index of list also starts with 0 and the list can be sliced, concatenated, etc.
#

list1

list1[-2]

list2[0:3]

list3[2]

#  - Unlike strings, it is possible to update a list by changing an element:

list1[2]

list1[2] = 2001 #replacing the item on position [2]  by the value 2001

list1[2]

list1

del list1[2] # deleting  item on position [2] from the list (2001)

list1

# # List operations
#
#  - The same as with strings for `+` (concatenation) and `*` (repetition) and most of operators. Ex: `len(list)`
#
#
#  -  More examples:
#

list1= ['physics', 'chemistry', 1997, 2000, 2004, 1999,2000]

len(list1)

2001 in list1


2000 in list1


list1.count(2000)

list1=[ str(s) for s in list1 ]

list1.sort()

list1

# # Hans-on 
#
# - From the list provided: 
#
#     1 - create a new list with the lengh of each item 
#
#     2 - Calculate the number of occurrences of the string "google"

list_hands_on = [11111, "python","perl","bash","R","google","microsoft","apple","google"]






# ---
# # Tuple
#
#  - Tuple is very similar to list, but immutable.
#
#
#  -  Instead of `[` `]` it uses `(` `)`
#

tup1 = ('bla','ble','blu',1,30)
tup2 = (1,2,3,4,5,6) 
tup3=(50)


# # Dictionary
#
#  - In the dictionary, for every item there is a `key` and a `value`, separated by  colon (`:`)
#  
#
#  - The items are separated by comma (`,`)
#  
#  
#  - Keys are unique in a dictionary, but the values may not be. 
#  
#
#  - You can access the values by the key 
#

dict1 = {'Name': 'John', 'Age': 25, 'Class': 'First'}

dict1['Name']


dict1['Age']

# # Comparison and logical operators in Python
#
# ## Logical:
#
# `and` -  both of both operators are true, the result is true
#
# `or`  - one of the operators is true, the result is true
#
# `not`  -  reverse the logical operator
#
# ## Comparison:
#
# `==`  check if two values are equal. Ex: `a == b`.
#
# `!=`  check if two values are NOT equal.
#
# `>`   check if the left is greater than the right. Ex: `a>b`. 
#
# `<`   check if the left is less than right
#
# `>=`  bigger or equal
#
# `<=`  smaller or equal
#

# # Decision making
# &nbsp;
#
# ![7](7.png)
#

#  - For decision making it is possible to use `if` just for one single condition, `if` and `elif` for more than one condition and `else` for anything else (and every elif). 
#

var1=100


if var1>90:
    print("Variable is greater than 90")

if var1>90:
    print("Variable is greater than 90")
else:
    print("Variable can be equal or smaller than 90")


var1=90

if var1>90:
    print("Variable is greater than 90")
else:
    print("Variable can be equal or smaller than 90")

if var1>90: 
    print("Variable is greater than 90")
elif var1<90:
    print("Variable is smaller than 90")
else:  #(could also be elif var1 == 90:)
    print("Variable is 90")

# # Identation
#
# ### Attention for the indentation in Python: 
#
#  - No braces or reserved words to indicate blocks of code for class, functions or flow control. It is done by indentation.
#
# ![8](8.png)
#

# # Loops
#
# - Sometimes you need to repeat the same operation, so it is possible to use a loop, with a structure of repetition based in a decision or just a fixed number of times.
#
# ![9](9.png)
#

#  -  `for` loop has the ability to iterate over the items of any sequence, such as a `list` or a `string`.
# ![10](10.png)
#

for letter in 'Python bla': 
    print('Current Letter :', letter)


fruits = ["apple", "mango", "banana"]
for fruit in fruits:
    print('Current fruit :', fruit)


fruits = ["apple", "mango", "banana"]
for fruit in fruits:
    for letter in fruit:
        print('Current fruit:', fruit, '\tcurrent letter :', letter)

#  - `while` loop statement in Python programming language repeatedly executes a target statement as long as a given condition is true.
#
# ![11](11.png)
#

count = 0 
while (count < 9):
    print("The count is:", count)
    count = count+1
print("Bye!")

# -  The infinite loop: you need to take care when using while loops to not create an endless loop  (unless you need one).
#
# ```
# var == 1 
# while var == 1:
#     print "Endless loop"
# ```
#
# -  To stop it you need to use CTRL+C
#

# # Hands on:
#
# - Create a list of words on the string "I'm learning python for bioinformatics"
# - Create a loop and print all words from the list, but just if there is "n" in the word







# # I/O in Python
#
#  - A program (or script) needs to deal with an input to generate an output. In Python there are many ways to input data and we will see the most basic ones first.
#
#

#  - The `input([prompt])` function reads one line from standard input and returns it as a string (removing the trailing newline).
#

str_1 = input("Enter your input: ")
print("Received input is : ", str_1)


#  - Before you read or write a file, you have to open it with the `open()` function.
#  - This function creates a file object:

afile = open("example.txt", "w+")


#  - The first argument is the file name, and the second is the “access mode”: It determines if the file will be only read, if can be write, append, etc. The default is r (read). 
#  - In the example, w+ is for both writing and reading, overwriting the existing file if the file exists. 
#
#  - Example of other modes are: 
#
# **‘w’** – Write mode which is used to edit and write new information to the file (any existing files with the same name will be erased when this mode is activated) 
#  
# **‘a’** – Appending mode, which is used to add new data to the end of the file; that is new information is automatically amended to the end 
#
# **‘r+’** – Special read and write mode, which is used to handle both actions when working with a file 
#
#

afile #afile is an object 

#  - Writing text on the new file

afile.write( "Python is amazing.\nYeah its great!!\n")  # \n is new line 


print(len("Python is amazing.\nYeah its great!!\n"))

#  - Closing the file:

afile.close()

# - Open, read and close a file:

fo = open("example.txt", "r+")


str_1 = fo.read()


type(str_1)

print(str_1.endswith)

str_1.splitlines()[1]

fo = open("example.txt", "r+")
str_2 = fo.read();
print(str_2)
fo.close()


fo.close()


# {{{
### include fasta file reading and manupulations 
# }}}





# # Functions
# &nbsp;
#
# - A function is a block of organized, reusable code that is used to perform a single, related action. Functions provide better modularity for your application and a high degree of code reusing.
#
#
# - As you already know, Python gives you many built-in functions like print(), etc. but you can also create your own functions. These functions are called user-defined functions.
#
#
# - You can define functions to provide the required functionality. 
#
#
# - Function blocks begin with the keyword `def` followed by the function name and parentheses `(` `)`:
#

def a_new_function():
    print("A new function")
    return


a_new_function()


# - Any input parameters or arguments should be placed within these parentheses. You can also define parameters inside these parentheses.
#
#
# - The first statement of a function can be an optional statement - the documentation string of the function or docstring.
#
#
# - The code block within every function starts with a colon (`:`) and is indented.
#
#
# - The statement `return [expression]` exits a function, optionally passing back an expression to the caller. A `return` statement with no arguments is the same as `return None`.
#

def a_new_function(a_str):
    print(a_str)
    return


a_new_function("another string")


# - It is possible to use arguments as keyword arguments. It allows you to skip arguments or place them out of order because python will identify by the keyword
#

# Function definition is here
def printinfo( name, age ):
    """
    This prints a passed info into this function
    """
    print("Name: ", name)
    print("Age ", age)
    return;


# Now you can call printinfo function - no keyword, but same order 
printinfo("miki",50)

printinfo(50,"miki")

# Now you can call printinfo function - using the keywords 
printinfo(age=50,name="miki")

# - How to get help:

help(printinfo)


# - It is possible to use default arguments:
#

def printinfo( name, age = 35 ): # default argument age is 35
    "This prints a passed info into this function"
    print("Name: ", name)
    print("Age ", age)
    return;



printinfo( name="miki" ) # no age argument, it will be the default


def sum_( arg1, arg2 ):
   # Add both the parameters and return them."
    total = arg1 + arg2
    print("Inside the function : ", total)
    return total


sum_( 10, 20 )


# ### Global vs. Local variables:
#
#  - Variables that are defined inside a function body have a local scope, and those defined outside have a global scope.
#

total = 0 # Global variable.


def sum_2( arg1, arg2 ):
    # Add both the parameters and return them."
    total = arg1 + arg2; # Local variable.
    print("Inside the function local total : ", total)
    return



sum_2( 10, 20 );
print("Outside the function global total : ", total )


# # Modules
# &nbsp;
#
#
# - A module allows you to logically organize your Python code. Grouping related code into a module makes the code easier to understand and use. 
#
#
# - A module is a Python object with arbitrarily named attributes that you can bind and reference.
#
#
# - Simply, a module is a file consisting of Python code. A module can define functions, classes and variables. A module can also include runnable code.
#
#
# - The Python code for a module named aname normally resides in a file named aname.py. Here's an example of a simple module, support.py
#

def print_func( name ):
    print("Hello : ", name)
    return



# {{{
# Creating a module on fly - just illustrative, never do it, use a text editor or IDE
# }}}

new_mod=open("support.py","w+")
new_mod.write("""\
def print_func( name ):\n\
\tprint(\"Hello : \", name)\n\
\treturn\
""")
new_mod.close()

# ### To use a module you first need to import it:

# {{{
# Import module support
import support

# Now you can call defined function that module as follows
support.print_func("Zara")

# }}}

# ### You can also use an alias for the module:

# Import module support with an alias “sup”
import support as sup


sup.print_func("Zara")


# # Hands-on
#
# - Write a Python function that takes a list of words and returns the length of the longest one.





# # Questions ???

# ![12](12.png)
#

# ---
# &nbsp;
# &nbsp;
#
#
# # Second Day workshop
# &nbsp;
#
#  - Numerical operations (Numpy)
#  
#  
#  - Dataframes (pandas)
#  
#  
#  - Plots (matplotlib and seaborn)
# &nbsp;
#
#

# # Numerical operations with Numpy
# &nbsp;
#
# - Numpy is the core library for scientific computing in Python. It provides a high-performance multidimensional array object, and tools for working with these arrays (including many numerical operations)
#
#
# - We will not cover arrays in detail, but we will see how to use numpy to do numerical operations
#
#
# - NumPy’s array class is called ndarray. It is also known by the alias array. Note that numpy.array is not the same as the Standard Python Library class array.array, which only handles one-dimensional arrays and offers less functionality.
#

import numpy as np

a = np.arange(15).reshape(3,5)

a

a = np.array([2,3,4])
b = np.array([2.1,3.2,4.5])

a

b

# # Operations with Numpy



c=a+b
c



d=a-b
d

d**2

a*b # Element wise product

a.dot(b) # Matrix product

e=np.arange(12).reshape(3,4)
e

e.sum(axis=0)

e.sum(axis=1)

e.min(axis=0)

e.min(axis=1)

e.max(axis=1)

e.cumsum(axis=1)


e.T

e[2,3]

e[-1]

# # Hands-on
#
# - with the 2 provided arrays:
#
# 1 - reshape both to 3 rows and 4 columns
#
# 2 - Transpose both
#
# 3 - calculate the sum of the cumsum of both arrays 
#
#

anewarray=np.random.rand(4,3)
anotherarray=np.random.rand(4,3)





# # Introduction to data frames with Pandas
#
# ### What is a Pandas data frame?
#
# - It is a two-dimensional size-mutable, potentially heterogeneous tabular data structure with labeled axes (rows and columns).
#
#
# - Arithmetic operations align on both row and column labels.
#
#

# ![13](13.png)
#

# - In general, we import a tabular or CSV file to a pandas dataframe
#
#
# - The are some options, for example, we can use different separator character, we can read a file with or without a header
#

# {{{
import pandas as pd


df=pd.read_csv("new_AVS.tsv")
# }}}

df.head()

#checking the first rows of dataframe - It shows 5 rows by default, but you can pass argument inside ()
df.head()

df.tail() # Last 5 rows 

#  - By default, function `read_csv()` would use `,` as separator to read the table. 
#
#
#  - We can change it for any carachter, for example tabular `\t`:

df=pd.read_csv("new_AVS.tsv",sep="\t")

df.head()

# - It's also possible to read excel tables using `read_excel("file.xlsx")`

df_cog=pd.read_excel("CAMERA_COG_RPKG_sig.xlsx")

# #### This table stores the normalized abundance of each COG (ortholog genes) on each sample

df_cog.head(10)

# # Working with the provided files
#
# ### 2 files from metagenomics analysis: 
#
# - `new_avs.csv` and `metatable.csv`
#
#
# - AVS stands for Average Genome Size calculated based on gene markers for each metagenomic sample
#
#
# - Metatable contains metadata about each sample
#

df=pd.read_csv("new_AVS.tsv",sep="\t") 

# ## This file is a parsed output from MicrobeCensus software:
#
# https://github.com/snayfach/MicrobeCensus
#
#
# MicrobeCensus is a fast and easy to use pipeline for estimating the average genome size (AGS) of a microbial community from metagenomic data.

df.head()  

# - The samples used on this study are from GOS - baltic sea. 

# ### Checking the column names:

df.columns

l=list(df.columns)

l

# ### Subseting a dataframe:

# ### We will start by looking at 2 columns only

new_df=df[["metagenome:","average_genome_size:","genome_equivalents:"]]

new_df.head()

# ### Rename columns (all): 

new_df.columns=["Sample","AGS","GE"]

new_df.head()

# ### Rename specific columns:

new_df.rename(columns= {"Sample":"sample"},inplace=False).head() #Needs inplace = True to really change it

# ### Working with index on dataframes:
#
#  - this will set the `index` and display the `new_df` with new `index`, without changing the original dataframe:

new_df.set_index("Sample").head() #Just show the dataframe with the index on "Sample", 
                                  # no alterations to original df

#Original still unchanged
new_df.head()

new_df.set_index("Sample",inplace=True) # inplace=True modifies the dataframe

# Now its changed
new_df.head()

# ### General statistic on dataframe:

new_df.describe()

# ### Get the min value from GE:

new_df["GE"].min()


# ### Get the min value from all columns:

new_df.min()


new_df.max()


new_df.dtypes

type(new_df)

# ### Transpose dataframe:

new_df.T

# {{{
#some hands on here
# }}}

# ### Sort by values in one column (in this case, GE):

new_df.sort_values(by="GE").head()


# ### By default it will sort values in an ascending fashin, but you can specify to be descending:

new_df.sort_values(by="GE",ascending=False).head()


# ### Selecting one column as Pandas Series:

new_df["AGS"].head()

# ### Selecting one column as new Pandas DataFrame:

new_df[["AGS"]].head()

# ### Slicing by rows (positional):

new_df.iloc[0:3]

# ### Slicing by rows (by index name):

new_df.loc["CAM_SMPL_003398":"CAM_SMPL_003418"]

# ### Conditional selection:

new_df[(new_df["GE"]<40) & (new_df["AGS"]>1000000)][["AGS"]].head()

new_df[new_df["GE"]>40] # Select the rows where GE is bigger than 40

# ### Operations with the values:

new_df.max()

new_df.min()

new_df.max()-new_df.min()

# ### Create a new column based on results of operations:

# {{{
# ind=new_df.index.tolist()
# new_df.reset_index(inplace=True, drop=True)
# new_df.index=ind
# new_df
# }}}

new_df.head()

new_df.loc[:,"GE - AGS"]=new_df.loc[:,"GE"]-new_df.loc[:,"AGS"]


new_df.head()

# ### Delete a column:

del new_df["GE - AGS"]

new_df.head()

new_df["AGS - GE"]=new_df["AGS"]-new_df["GE"]


# #### More about Pandas in: http://pandas.pydata.org/
#
# #### Cheat Sheet: http://pandas.pydata.org/Pandas_Cheat_Sheet.pdf

# # Hands on:
#
#
# - read a dataframe (any!! From internet, from your computer...)
#
#
# -  print the columns 
#
#
# -  store the column names into a list
#
#
# - print the max value from the last column 
#
#
# - print the max value from the last row 
#
#
#





# ---
# # Plots with matplotlib

# - Matplotlib is a Python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms. Matplotlib can be used in Python scripts, the Python and IPython shells, the Jupyter notebook, web application servers, and four graphical user interface toolkits.
#
# #### https://matplotlib.org/

import matplotlib.pyplot as plt

# - Seaborn is a Python visualization library based on matplotlib. It provides a high-level interface for drawing attractive statistical graphics.
#
# #### https://seaborn.pydata.org/

import seaborn as sns

# ### Starting with plots - the ugly and lazy way

#function plot direct from pandas dataframe
new_df["GE"].plot()
plt.show()

# {{{
#Rotate the name of samples 

new_df["GE"].plot()


plt.xticks(rotation=90)

plt.show()

# }}}

# ### This plot is useless, we dont have continuos data, we need a bar plot in this case

new_df["GE"].plot(kind="bar")  # kind can be also pie, scatter, etc ...
plt.xticks(rotation=90)
plt.show()

# ### We need to control the size  and style of figure, and also axis titles, etc

# {{{
sns.set_style("white")

#Define a figure and its size
plt.figure(figsize=(15,8))



#Bar plot 

# range(len(new_df["AGS"])) -> define the number of bars to be plotted and the position on X axis 
x=range(len(new_df["AGS"]))

# new_df["AGS"] -> values for axis Y
y=new_df["AGS"]

plt.bar(x,y,color="b")

#axis labels
plt.ylabel("AGS ",fontsize="16")
plt.xlabel("Samples",fontsize="16")

plt.title("PLOT FOR WORKSHOP",fontsize="22")
#xticks -> define the text to be ploted on axis X -> here we use index (the name of samples)

# plt.xticks(range(0,len(new_df.index),2), new_df.index,rotation="vertical",fontsize="12")
### error .... index values not matching the positions of xticks

plt.xticks(range(0,len(new_df.index),2), list(new_df.index)[0::2],rotation="vertical",fontsize="12")
plt.yticks(fontsize="22")

## Save the figure in pdf

plt.savefig("my_first_plot.pdf", format = 'pdf', dpi = 300, bbox_inches = 'tight')

#Necessary to display the plot here
plt.show()
# }}}

# #### More about seaborn styles in: https://seaborn.pydata.org/generated/seaborn.set_style.html

# ### The values on Y are in bases, hard to read, so we can convert to Mb

new_df["AGS (mb)"]=new_df["AGS"]/1000000

# {{{
#Seaborn style
sns.set_style("white")

#Define a figure and its size
plt.figure(figsize=(15,8))


x=range(len(new_df["AGS (mb)"]))

y=new_df["AGS (mb)"]


#Bar plot 
plt.bar(x,y,color="b")

#axis labels
plt.ylabel("AGS (Mb)",fontsize="16")
plt.xlabel("Samples",fontsize="16")


plt.xticks(range(len(new_df.index)), new_df.index,rotation="vertical")

plt.show()
# }}}

# ### Metadata from samples

df2=pd.read_csv("metatable.csv", encoding = "latin")


df2.head()

new_df.head()

new_df.reset_index(inplace=True)

new_df=new_df.rename(columns={"index":"Sample"})

new_df.head()

# {{{
# checking size of both dataframes, if they have same number of samples
# }}}

len(new_df)

len(df2)

#mergeing 2 df
merged=pd.merge(new_df,df2,left_on="Sample",right_on="SAMPLE_ACC")


merged.head(15)

l1=list(merged.columns)

#removing columns - (axis 1) -  where ALL the values are NaN (null, missing)
merged.dropna(how="all",axis=1,inplace=True)


merged.head()

l2=list(merged.columns)

[x for x in l1 if x not in l2]


# ### Sort values by filter size:

merged=merged.sort_values("FILTER_MIN")

merged.head(10)

merged.reset_index(inplace=True,drop=True)

merged.head(10)

# ### Creating 3 groups (g1, g2, g3) for each filter_min size

merged["FILTER_MIN"].drop_duplicates()

g1_ = merged[merged["FILTER_MIN"]==0.1]
g2_ = merged[merged["FILTER_MIN"]==0.8]
g3_ = merged[merged["FILTER_MIN"]==3.0]


# {{{
# Group by and calculating the mean of each group
# it will correspond to the same order we saw on drop_duplicates

list(merged.groupby("FILTER_MIN")["AGS (mb)"].mean())
# }}}

# # Hands-on:
#
# - plot a barplot with GE columns only for g1 (FILTER_MIN = 0.1).
# - Add the legend of axis x (xticks) with 45 degrees rotation, using name of samples





# ### Creating more complicated bar plots:
#
# - ploting 3 groups (g1, g2 and g3), different colors, and a line defining the average value for each group

# {{{
sns.set_style("white")

plt.figure(figsize=(15,10))

#Here we are ploting on the same figure, 3 bar plots, ranging from the original index obtained from the merged dataframe

plt.bar(g1_.index,g1_["AGS (mb)"],color="b") # -> color b=blue, each bar plot we use different colors

plt.bar(g2_.index,g2_["AGS (mb)"],color="c")

plt.bar(g3_.index,g3_["AGS (mb)"],color="r")

plt.ylabel("AGS (mb)",size="16") 
plt.xticks(range(len(merged)), merged["Sample"], size='small',rotation="vertical")
plt.legend(merged["FILTER_MIN"].drop_duplicates())

#here we plot 3 horizontal lines (hlines) with the mean of AGS values for each group
plt.hlines(g1_["AGS (mb)"].astype(float).mean(),g1_.index[0],g1_.index[-1])
plt.hlines(g2_["AGS (mb)"].astype(float).mean(),g2_.index[0],g2_.index[-1])
plt.hlines(g3_["AGS (mb)"].astype(float).mean(),g3_.index[0],g3_.index[-1])
plt.show()

# }}}

# ### Starting with seaborn, using hue and crating more fancy plots

# {{{
sns.set_style("white")

plt.figure(figsize=(15,10))
x=merged.index
y="AGS (mb)"
sns.barplot(x=x,y=y,hue="FILTER_MIN",data=merged)
plt.xticks(range(len(x)), merged["Sample"], size='small',rotation="vertical")

plt.show()
# }}}

# # Hands-on: 
#
# - add the mean values line to the previous plot





# ### Barplot with error bars, violinplot and boxplots

# {{{
sns.set_style("white")

plt.figure(figsize=(15,10))
x=merged["FILTER_MIN"]
y=merged["AGS (mb)"]
sns.barplot(x=x,y=y)
plt.show()
# }}}

# {{{
sns.set_style("white")

plt.figure(figsize=(15,10))
x=merged["FILTER_MIN"]
y=merged["AGS (mb)"]
sns.violinplot(x=x,y=y)
plt.show()
# }}}

# {{{
sns.set_style("white")

plt.figure(figsize=(15,10))
x=merged["FILTER_MIN"]
y=merged["AGS (mb)"]
sns.boxplot(x=x,y=y)
plt.show()
# }}}



# # Hands on:
#
# - create the same violin plots, but adding dots on the plot
# - create violin plots for sample volume, instead filter min





# ## Exploring data, checking correlations

# {{{
#correlation between 2 variables

merged[["AGS (mb)","temperature - (C)"]].corr()
# }}}

help(merged[["AGS (mb)","temperature - (C)"]].corr)

# ## Do you think that AGS (average genome size) correlates with temperature?

# {{{
#Doing a scatter plot

sns.set_style("white")
x=merged["AGS (mb)"]
y=merged["temperature - (C)"]
plt.figure(figsize=(15,10))
plt.scatter(x,y)
plt.ylabel("temperature - (C)",size="16")
plt.xlabel("AGS (mb)",size="16")
plt.show()
# }}}

sns.jointplot(x="temperature - (C)", y="AGS (mb)", data=merged, kind="reg",
                  xlim=(4, 20), ylim=(0.5, 2.5), color="r")
plt.show()

sns.jointplot(x="LATITUDE", y="AGS (mb)", data=merged, kind="reg",
                  xlim=(10, 20), ylim=(0.5, 3), color="r")
plt.show()

# # Hands on: 
# - plot correlation plot with other variables





# ## Heatmaps
#
# - Working with genes abundance on samples

df_cog.index=df_cog["Unnamed: 0"].tolist()
df_cog=df_cog.drop(["Unnamed: 0"], axis=1)

# Remove rows only with 0
df_2=df_cog[(df_cog.T != 0).any()]

df_2.columns

len(df_cog) # number of cogs (gene ortologs groups) on the original table

len(df_2) # number of cogs after removing those with 0 abundance in all the samples

# {{{
# New dataframe only with the RPKG for all the samples

heat=df_2[[x for x in df_2.columns if "CAM_SMPL"in x ]]
# }}}

heat.head()

# {{{
##test.ix[:, test.columns != 'compliance']
# }}}

sns.clustermap(heat,figsize=(12,25),cbar_kws={'label': 'COGs RPKG'},yticklabels=False)
#plt.savefig("../clustermap_cogs.pdf", format = 'pdf', dpi = 300, bbox_inches = 'tight')
plt.show()

# {{{
# customizing colors 
# }}}

import matplotlib.cm as cm
from  matplotlib.colors import LinearSegmentedColormap

c = ["pink","lightcoral","red","darkred","black"]
v = [0,.25,.5,.75,1.]
l = list(zip(v,c))
cmap=LinearSegmentedColormap.from_list('rg',l, N=1024)

sns.clustermap(heat,figsize=(12,25),cbar_kws={'label': 'COGs RPKG'},yticklabels=False,cmap=cmap)
#plt.savefig("../clustermap_cogs.pdf", format = 'pdf', dpi = 300, bbox_inches = 'tight')
plt.show()

heat_log= heat.apply(lambda x: np.log10(x+0.001))

c = ["darkred","red","lightcoral","white","palegreen","green","darkgreen"]
v = [0,.15,.4,.5,0.6,.9,1.]
l = list(zip(v,c))
cmap=LinearSegmentedColormap.from_list('rg',l, N=1024)

sns.clustermap(heat_log,figsize=(12,25),cbar_kws={'label': 'log10(COGs RPKG)'},yticklabels=False,cmap=cmap)
#plt.savefig("../clustermap_cogs.pdf", format = 'pdf', dpi = 300, bbox_inches = 'tight')
plt.show()

heat_log= heat.apply(lambda x: np.log10(x+0.001))
vals=heat_log.values
vals=list(vals)
vals=list([ list(v) for v in vals ] )
vals=[ item for row in vals for item in row]
kde=sns.kdeplot(vals)

# {{{
x=list(list(kde.get_lines()[0].get_data())[0])
y=list(list(kde.get_lines()[0].get_data())[1])
xy=pd.DataFrame({"x":x,"y":y})
# #xy.head()
miny=max( xy[ xy["x"]<-2 ]["y"].tolist() )
maxy=max( xy[ xy["x"]>-2 ]["y"].tolist() )
minx=xy[xy["y"]==miny]["x"].tolist()[0]
maxx=xy[xy["y"]==maxy]["x"].tolist()[0]

centery=min( xy[ ( xy["x"]>minx )  & ( xy["x"]<maxx ) ]["y"].tolist() )
centerx=xy[xy["y"]==centery]["x"].tolist()[0]
print(minx,centerx,maxx)
# }}}

# {{{
sns.clustermap(heat_log,figsize=(12,25),cbar_kws={'label': 'log10(COGs RPKG)'},yticklabels=False,cmap='bwr', vmin=minx, center=centerx ,vmax=maxx)
# sns.clustermap(heat_log,figsize=(12,25),cbar_kws={'label': 'log2(COGs RPKG)'},yticklabels=False,cmap='bwr', vmin=minx, center=0 ,vmax=maxx)

#plt.savefig("../clustermap_cogs.bwr.pdf", format = 'pdf', dpi = 300, bbox_inches = 'tight')
plt.show()
# }}}

sns.clustermap

# # Hands-on:
#
# - Plot the same heatmaps, but only for genes where pvalue adjusted < 0.01
# - plot a heatmap for z-score 





# ## Getting upregulated COGs  in the comparisons of different filtrations methods groups

a=set(df_2[df_2["log2(FC)_g1/g3"]>0].index)
b=set(df_2[df_2["log2(FC)_g1/g2"]>0].index)
c=set(df_2[df_2["log2(FC)_g2/g3"]>0].index)

# ### To install matplotlib_venn:
#
# ```pip install --user matplotlib-venn```
#
# - more info about the package in https://pypi.org/project/matplotlib-venn/

from matplotlib_venn import venn3,venn3_unweighted 


# {{{
plt.figure(figsize=(10,10))


vd = venn3_unweighted([a, b,c], ('up g1/g3', 'up g1/g2','up g2/g3'))

plt.show()
# }}}

# {{{
plt.figure(figsize=(10,10))


vd = venn3([a, b,c], ('up g1/g3', 'up g1/g2','up g2/g3'))

plt.show()
# }}}

# # hypergeometric test
#
# The hypergeometric test uses the hypergeometric distribution to measure the statistical significance of having 
# drawn a sample consisting of a specific number of k successes (out of n total draws) from a population of size 
# N containing K successes. In a test for over-representation of successes in the sample, the hypergeometric 
# p-value is calculated as the probability of randomly drawing k or more successes from the population in n total draws. 
# In a test for under-representation, the p-value is the probability of randomly drawing k or fewer successes.
#
# Biologist and statistician Ronald Fisher
# The test based on the hypergeometric distribution (hypergeometric test) is identical to the corresponding 
# one-tailed version of Fisher's exact test.[6] Reciprocally, the p-value of a two-sided Fisher's exact test 
# can be calculated as the sum of two appropriate hypergeometric tests (for more information see[7]).
#
# https://en.wikipedia.org/wiki/Hypergeometric_distribution

# {{{
from scipy.stats import hypergeom

M=22000 # total number of genes in organism / population size (previously N)
n=5000 # genes in group I / number of successes in the population (previously K)
N=3000 # genes in group II / sample size (previously n)
x=1000 # intersect / number of drawn “successes” (previously k)

p=hypergeom.sf(x-1, M,n,N)
print(p)
# }}}

# # Hands-on
#
# - Plot a heatmap for log2(FC) up regulated genes for the 3 comparions (g1/g2, g2/g3 and g1/g3)





# ### Package sklearn - more info at  http://scikit-learn.org/stable/index.html

# # PCA - Principal component analysis
#
# Principal component analysis (PCA) is a linear dimensionality reduction technique with applications in exploratory data analysis, 
# visualization and data preprocessing. The data is linearly transformed onto a new coordinate system such that the directions 
# (principal components) capturing the largest variation in the data can be easily identified.
#
# The principal components of a collection of points in a real coordinate space are a sequence of 
# p unit vectors, where the i-th vector is the direction of a line that best fits the data while being orthogonal to the first 
# i−1 vectors. Here, a best-fitting line is defined as one that minimizes the average squared perpendicular distance from the points 
# to the line. These directions (i.e., principal components) constitute an orthonormal basis in which different individual dimensions 
# of the data are linearly uncorrelated. Many studies use the first two principal components in order to plot the data in two 
# dimensions and to visually identify clusters of closely related data points.[1]
#
# https://en.wikipedia.org/wiki/Principal_component_analysis
#
# Very nice expanation about PCA: https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues
#
# [hbctraining.github.io](https://hbctraining.github.io/DGE_workshop/lessons/principal_component_analysis.html)

from sklearn import preprocessing
from sklearn.decomposition import PCA
from itertools import cycle

heat.head(2)

len(heat)

# {{{
# we have 164 COGs, that means, we have 164 variables to be transformed to 2 (or more) dimensions
# }}}

df_pca=heat.T.reset_index()

df_pca.head()

merged.head()

df_pca=pd.merge(merged[["Sample","FILTER_MIN"]],df_pca,left_on="Sample",right_on="index").drop("index",axis=1)

df_pca.set_index(["Sample","FILTER_MIN"],inplace=True)


df_pca.head()

# ### feature scaling
#
# [scikit-learn.org](https://scikit-learn.org/stable/auto_examples/preprocessing/plot_scaling_importance.html)
#
# Feature scaling through standardization, also called Z-score normalization, is an important preprocessing step for many machine learning algorithms. It involves rescaling each feature such that it has a standard deviation of 1 and a mean of 0.
#
# Even if tree based models are (almost) not affected by scaling, many other algorithms require features to be normalized, often for different reasons: to ease the convergence (such as a non-penalized logistic regression), to create a completely different model fit compared to the fit with unscaled data (such as KNeighbors models). 

# {{{
pca = PCA(copy=True, iterated_power='auto', n_components=2, random_state=None,
  svd_solver='auto', tol=0.0, whiten=False)

#scaling the values
df_pca_scaled = preprocessing.scale(df_pca)

projected=pca.fit_transform(df_pca_scaled)
# }}}

print(pca.explained_variance_ratio_)


tmp=pd.DataFrame(projected)
tmp.rename(columns={0: 'Component 1',1: 'Component 2'}, inplace=True)
tmp.head()

# {{{
final_pca=pd.merge(df_pca.reset_index(),tmp,left_index=True,right_index=True)

# in the final_pca dataframe, in adition to all the COGs, 
# we have new columns (0 and 1), with the PCA results (2 components)
final_pca.head()
# }}}

# {{{
sns.set_style("white")
color_gen = cycle(('blue', 'green', 'red'))
plt.figure(figsize=(10,10))
for lab in set(final_pca["FILTER_MIN"]):
    plt.scatter(final_pca.loc[final_pca['FILTER_MIN'] == lab, 'Component 1'], 
                final_pca.loc[final_pca['FILTER_MIN'] == lab, 'Component 2'], 
                c=next(color_gen),
                label=lab)

plt.xlabel('component 1  - ' + str(pca.explained_variance_ratio_[0]*100)+ " % " )
plt.ylabel('component 2  - '+ str(pca.explained_variance_ratio_[1]*100)+ "  % " )
plt.legend(loc='best')
plt.show()
# }}}

# {{{
#pca.components_
# }}}

# Principal axes in feature space, representing the directions of maximum variance in the data. 
components = pd.DataFrame(pca.components_, columns = df_pca.columns, index=[1, 2])
components


# ## Visualize Loadings
#
# It is also possible to visualize loadings using shapes, and use annotations to indicate which feature a certain 
# loading original belong to. Here, we define loadings as:
#     
# ```
# loadings=eigenvectors * sqrt(eigenvalues)
# ```
#

loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
loadings = pd.DataFrame(loadings, index=df_pca.columns.values, columns=["c1","c2"] )
for c in loadings.columns.tolist():
    vals=loadings[c].tolist()
    vals=[min(vals),max(vals)]
    loadings.loc[ loadings[c].isin(vals) , f"key {c}"] = "yes"
important_vectors=loadings[ (loadings["key c1"]=="yes") | (loadings["key c2"]=="yes") ]
important_vectors

# {{{
color_gen = cycle(('red', 'green', 'blue'))
plt.figure(figsize=(10,10))
for lab in set(final_pca["FILTER_MIN"]):
    plt.scatter(final_pca.loc[final_pca['FILTER_MIN'] == lab, 'Component 1'], 
                final_pca.loc[final_pca['FILTER_MIN'] == lab, 'Component 2'], 
                c=next(color_gen),
                label=lab,alpha=0.8)

vectors=loadings[ (loadings["key c1"]=="yes") | (loadings["key c2"]=="yes") ]
for i in list(vectors.index):
    x=vectors.loc[i,"c1"]
    y=vectors.loc[i,"c2"]
    # s=vectors.loc[i,"c1xc2"]
    plt.arrow(0, 0, x, y, color='k', width=0.0005, head_width=0.25, alpha=0.75)
    plt.text(x*1.2, y*1.2, i, color='k', alpha=0.75) 

plt.xlabel('Component 1  - ' + str(pca.explained_variance_ratio_[0]*100)+ " % " )
plt.ylabel('Component 2  - '+ str(pca.explained_variance_ratio_[1]*100)+ "  % " )
plt.legend(loc='best')
plt.show()
# }}}

# # Hands-on:
#
# - Plot the same PCA, but fixing the text on arrows (keeping only COG number, removing the descriptions)
# - Add a vertical line and a horizontal line on zero values 







# ### Plotting the 4 COGs with extreme vectors on PCA

cols=list(important_vectors.index)

sns.set_style("white")
plt.figure(figsize=(15,10))
tmp=final_pca[final_pca['Sample'].isin(g1_["Sample"])]
tmp2=final_pca[final_pca['Sample'].isin(g2_["Sample"])]
tmp3=final_pca[final_pca['Sample'].isin(g3_["Sample"])]
for c in cols:
    plt.figure(figsize=(15,10))
    plt.bar(tmp.index,tmp[c],color="b")
    plt.bar(tmp2.index,tmp2[c],color="g")
    plt.bar(tmp3.index,tmp3[c],color="r")
    plt.ylabel(c,size="12") 
    plt.xticks(range(len(final_pca)), final_pca["Sample"], size='small',rotation="vertical")
    plt.legend(final_pca["FILTER_MIN"].drop_duplicates())

    #here we plot 3 horizontal lines (hlines) with the mean of AGS values for each group
    plt.hlines(tmp[c].astype(float).mean(),g1_.index[0],g1_.index[-1])
    plt.hlines(tmp2[c].astype(float).mean(),g2_.index[0],g2_.index[-1])
    plt.hlines(tmp3[c].astype(float).mean(),g3_.index[0],g3_.index[-1])
    plt.show()


# {{{
sns.set_style("white")

for c in cols:
    plt.figure(figsize=(15,10))
    x=final_pca["FILTER_MIN"]
    y=final_pca[c]
    sns.violinplot(x=x,y=y)
plt.show()
# }}}

# {{{
## Kmeans

## clusterization (kmeans with euc distance)
## calculate distance btw centroid of clusters and each gene
## plot bar plots for the genes with shorter and bigger distance from one cluster

### a bit more about K-means: https://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/
# }}}

def dist(a, b, ax=1):
    return np.linalg.norm(a - b, axis=ax)


f1 = final_pca["Component 1"].values
f2 = final_pca["Component 2"].values
X = np.array(list(zip(f1, f2)))
plt.scatter(f1, f2, c='black', s=15)
plt.show()

from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist

# {{{
# k means determine k
distortions = []
K = range(1,10)
for k in K:
    kmeanModel = KMeans(n_clusters=k).fit(X)
    kmeanModel.fit(X)
    distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])

# Plot the elbow
plt.plot(K, distortions, 'bx-')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('The Elbow Method showing the optimal k')
plt.show()
# }}}

# {{{
from sklearn.cluster import KMeans


## scalled 

# Number of clusters
k=3
kmeans_ = KMeans(n_clusters=k)

# scale
X_=preprocessing.scale(X)

# Fitting the input data
kmeans_ = kmeans_.fit(X_)

# Getting the cluster labels
labels_ = kmeans_.predict(X_)

# Centroid values
C_ = kmeans_.cluster_centers_
print(C_)



## non scalled

# Number of clusters
k=3
kmeans = KMeans(n_clusters=k)

# # scale
# X_=preprocessing.scale(X)

# Fitting the input data
kmeans = kmeans.fit(X)

# Getting the cluster labels
labels = kmeans.predict(X)

# Centroid values
C = kmeans.cluster_centers_
print(C)

# }}}

from copy import deepcopy

# {{{

# colors = ['#FF0000', 'g', '#FF7D40', '#66CDAA', '#473C8B', '#9400D3','#8B8386','#FFC0CB','#B0171F']
# fig, ax = plt.subplots()
# for i in range(k):
#         points = np.array([X[j] for j in range(len(X)) if clusters[j] == i])
#         ax.scatter(points[:, 0], points[:, 1], s=15, c=colors[i])
# ax.scatter(C[:, 0], C[:, 1], marker='*', s=100, c='#050505')
# plt.show()
# }}}
for i in range(len(labels)):
    print(labels_[i], labels[i])

# {{{
color_gen = cycle(('red', 'green', 'blue'))
final_pca["KMeans (scaled)"]=labels_
final_pca["KMeans (non scaled)"]=labels

# plt.figure()

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5),sharey=True, sharex=True)
# fig.suptitle('Horizontally stacked subplots')
# ax1.plot(x, y)
# ax2.plot(x, -y)

color_gen = cycle(('red', 'green', 'blue'))
for lab in set(final_pca["FILTER_MIN"]):
    ax1.scatter(final_pca.loc[final_pca['FILTER_MIN'] == lab, 'Component 1'], 
                final_pca.loc[final_pca['FILTER_MIN'] == lab, 'Component 2'], 
                c=next(color_gen),
                label=lab,alpha=0.8)
ax1.set_title('FILTER_MIN')


color_gen = cycle(('red', 'green', 'blue'))
for lab in set(final_pca['KMeans (scaled)']):
    ax2.scatter(final_pca.loc[final_pca['KMeans (scaled)'] == lab, "Component 1"], 
                final_pca.loc[final_pca['KMeans (scaled)'] == lab, "Component 2"], 
                c=next(color_gen),
                label=lab,alpha=0.8)
ax2.set_title('KMeans (scalled)')

color_gen = cycle(('red', 'green', 'blue'))
for lab in set(final_pca['KMeans (non scaled)']):
    ax3.scatter(final_pca.loc[final_pca['KMeans (non scaled)'] == lab, "Component 1"], 
                final_pca.loc[final_pca['KMeans (non scaled)'] == lab, "Component 2"], 
                c=next(color_gen),
                label=lab,alpha=0.8)
ax3.set_title('KMeans (non scaled)')

for ax in [ax1,ax2,ax3]:
    ax.set(xlabel=f'PC1 {round(pca.explained_variance_ratio_[0]*100, 2)}%', 
           ylabel="")
    o=[x.set_linewidth(1.5) for x in ax.spines.values()]
    ax.legend()
    
ax1.set(xlabel=f'PC1 {round(pca.explained_variance_ratio_[0]*100, 2)}%', 
       ylabel=f'PC2 {round(pca.explained_variance_ratio_[1]*100, 2)}%')

plt.ylim(-18, 18)
plt.xlim(-18, 18)
plt.show()
# }}}

from sklearn.metrics import pairwise_distances_argmin_min


closest_, _ = pairwise_distances_argmin_min(kmeans_.cluster_centers_, X_)
print(list(closest_))

final_pca.loc[list(closest_)]


# ## Hands-on:
#
# - Try the kmeans clusterization with different number of cluters (4 and 5)

#
# # Lifespan analysis 

# - lifelines package: http://lifelines.readthedocs.io/en/latest/
# - [databricks tutorial](https://notebooks.databricks.com/notebooks/CME/Survival_Analysis/index.html#Survival_Analysis_1.html)

# {{{
from lifelines.datasets import load_waltons
df = load_waltons() # returns a Pandas DataFrame

print(df.head())


T = df['T']
E = df['E']
# }}}

# - T is an array of durations, E is a either boolean or binary array representing whether the “death” was observed (alternatively an individual can be censored).

print(list(E))

from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()
kmf.fit(T, event_observed=E)  # or, more succiently, kmf.fit(T, E)

kmf.survival_function_


kmf.median_survival_time_

kmf.plot(show_censors=True)
plt.show()

# {{{
groups = df['group']
ix = (groups == 'miR-137')

kmf.fit(T[~ix], E[~ix], label='control')
ax = kmf.plot()

kmf.fit(T[ix], E[ix], label='miR-137')
kmf.plot(ax=ax)
plt.show()
# }}}

# # log-rank test
#
# [https://www.ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5448258/#:~:text=The%20log%2Drank%20test%20is%20a%20nonparametric%20hypothesis%20test%20to,is%20time%20to%20an%20event.)
#
# The log-rank test is a nonparametric hypothesis test to compare the survival trend of two or more groups when there are censored observations. It is widely used in clinical trials to compare the effectiveness of interventions when the outcome is time to an event. 
#
# The null hypothesis for the test is that that there is no difference in the survival experience of the subjects in the different groups being compared. Its name derives from its relation to a test that uses the logarithms of the ranks of the data. There are certain preconditions for applying the test.
#
# - The test assumes no particular distribution for the survival curve, that is, it is distribution free (nonparametric)
# - Subjects who are censored have the same probability of the event as subjects who are fully followed up, that is, the censoring must be noninformative
# - The proportional hazards assumption must be met, that is, there is no tendency for one group to have better survival than the other group at earlier time points and then worse survival at later time points. Such a tendency would be reflected in Kaplan–Meier plots that diverge initially and then cross.
#
# [https://lifelines.readthedocs.io](https://lifelines.readthedocs.io/en/latest/lifelines.statistics.html#lifelines.statistics.logrank_test)
#
# Measures and reports on whether two intensity processes are different. That is, given two event series, determines whether the data generating processes are statistically different. The test-statistic is chi-squared under the null hypothesis
#
# Also: https://www.databricks.com/notebooks/survival_analysis/survival_analysis_02_exploratory_analysis.html

# {{{
from lifelines.statistics import logrank_test

results = logrank_test(T[ix], T[~ix], E[ix], E[~ix], alpha=.99)

results.print_summary()
# }}}

# {{{
## https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html

from lifelines.datasets import load_regression_dataset
regression_dataset = load_regression_dataset()

regression_dataset.head()
# }}}

# # Univariable Cox Regression 
#
# ## log-rank vs cox regression
#
# - log-rank test does not work for continuous exposures
# - it does not allow for covariate adjustment
# - the usual P-value from log-rank may not be as accurate as the likelihood ratio χ2 statistic from Cox PH
#
#

# {{{
from lifelines import CoxPHFitter

# Using Cox Proportional Hazards model
cph = CoxPHFitter()
cph.fit(regression_dataset, 'T', event_col='E',robust=True)
cph.print_summary()
# }}}

# ****
# - partial log-likelihood: The partial log-likelihood for a Cox regression model is simply the logarithm of the partial likelihood. In general, when comparing two models fit to the same data, the model with the larger log likelihood is considered a better “fit”. Note that these log likelihood values are often negative! In this case, the larger value is the same as the less negative value. Thus the model with the less negative value is considered the model with a better "fit".
# ****
# - coef: eg. 0.22 == 22% more
# - exp(coef) : Hazard Ratios; risk change by factor of x (1.25x)
# - cmp to
# - z: test statistic of the z-test ( Z-value is a test statistic for Z-tests that measures the difference between an observed statistic and its hypothesized population parameter in units of the standard deviation)
# - p: p-value
# - log2(p)
# ****
# - concordance: with Cox Assumptions
# - Partial AIC: This value is derived from an information theory approach that attempts to determine how well the data fit the model. A model with a small AIC value suggests a better fit than a model with a large AIC value on the same data. The last very important concept to note about AIC values is that they can only be compared between models fit to the same data! AIC values are calculated from likelihood, which is specific to the data set being analyzed.
# - log-likelihood ratio test	: 𝑋2  statistic from a likelihood ratio test; "Good" models result in a higher value of likelihood while "poor" models have a lower value for their likelihood.
# ****

np.exp(0.22)

cph.plot()
plt.show()

# # system report

# {{{
import pkg_resources
from datetime import datetime
import sys

# list installed packages
installed_packages = pkg_resources.working_set
installed_packages_list = sorted(["%s==%s" % (i.key, i.version)
   for i in installed_packages])

# get python version
v=sys.version

# write to file with time stamp
now = datetime.now()
dt_string = now.strftime("%d.%m.%Y-%H.%M.%S")
with open(f"packages.{dt_string}.txt", "w") as report:
    report.write(f"Python: {v}\n")
    report.write("\n".join(installed_packages_list))

# print versions
print(v)  
output=[ print(p) for p in installed_packages_list ]

# }}}
# ---
# &nbsp;
# &nbsp;
#
#
# # Third Day workshop
#
#  - biomart and gene annotation tables
#  - querying DAVID with Python
#  - hands-on large data: exploring GTEX RNAseq data, finding gene-gene correlations and age related changes
#     - Anova
#     - linear models (simple and multiple)
#  - image analysis: identifying and quantifying objects in microscopy pictures
#
#


# ## gene annotation tables
#
# Gene annotations tables can be downloaded from 
#
# https://www.ensembl.org/info/data/ftp/index.html
#
# We will download the gene annotation for Caenorhabditis elegans:
#
# https://ftp.ensembl.org/pub/release-111/gtf/caenorhabditis_elegans

import urllib
import gzip
url = 'https://ftp.ensembl.org/pub/release-111/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.111.gtf.gz'
filehandle, _ = urllib.request.urlretrieve(url)
with gzip.open(filehandle, 'rb') as f:
    file_content = f.read()
print(type(file_content))
file_content=file_content.decode("utf-8")
print(type(file_content))

file_content[:1000]

print(file_content[:1000])

gtf=file_content.split("\n")
gtf=[ s for s in gtf if s != "" ]
gtf=[ s.split("\t") for s in gtf if s[0] != "#" ]
gtf[:3]

gtf=pd.DataFrame(gtf,columns=['seqname','source','feature','start','end','score','strand','frame','attribute'])
gtf.head()



# {{{
## get table of gene names and gene ids

name_ids=gtf[gtf["feature"]=="gene"]
name_ids=name_ids.reset_index(inplace=False, drop=False)

def get_attribute(a, attributes):
    l=attributes.split(";")
    l=[ s.split(" ") for s in l]
    res=np.nan
    for s in l:
        if a in s:
            if '"' in s[-1]:
                res=s[-1][1:-1]
            else:
                res=s[-1]
    return res


name_ids["gene_name"]=name_ids["attribute"].apply(lambda x: get_attribute("gene_name", x) )
name_ids["gene_id"]=name_ids["attribute"].apply(lambda x: get_attribute("gene_id", x) )
name_ids     
# }}}

name_ids.loc[46921,"attribute"]

name_ids=name_ids[["gene_name","gene_id"]]
name_ids.head()

# # biomart

from biomart import BiomartServer
from io import StringIO
import pandas as pd

help(BiomartServer)

server=BiomartServer("http://www.ensembl.org/biomart")
help(server)

datasets=server.show_datasets()

datasets

help(server.show_datasets)


# {{{
# solution
# write a class that captures stdout
# and parses specifically biomart stdout
# A Class is like an object constructor, or a "blueprint" for creating objects.
# more on objects: https://www.w3schools.com/python/python_classes.asp

class BiomartOutput(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(  [ s.lstrip().rstrip(",").lstrip("{").rstrip("}") for s in self._stringio.getvalue().replace("'","").splitlines() ] )
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

with BiomartOutput() as datasets:
    server.show_datasets()

print("1st 10 datasets:\n", datasets[:10], "\n")
print("'sapiens' datasets:\n", "\n".join( [ s for s in datasets if "sapiens" in s ] ), "\n" )
print("'elegans' datasets:\n", "\n".join( [ s for s in datasets if "elegans" in s ] ), "\n" )
# }}}

# {{{
dataset="celegans_gene_ensembl"
organism=server.datasets[dataset]

with BiomartOutput() as attributes :
    organism.show_attributes()

attributes
# }}}

attributes_filtered=[ s for s in attributes if "homolog" not in s ]
attributes_filtered

query_attributes=[ "ensembl_gene_id", "external_gene_name", "gene_biotype",  "go_id", "name_1006"]
response=organism.search({"attributes":query_attributes})
bmdf=response.content.decode().split("\n")
bmdf=[s.split("\t") for s in bmdf ]
bmdf=pd.DataFrame(bmdf, columns=query_attributes)
bmdf

bmdf=bmdf.groupby(["ensembl_gene_id","external_gene_name","gene_biotype"], as_index=False).agg( {"go_id": "; ".join, "name_1006": "; ".join} )
bmdf.head()

bmdf[bmdf["external_gene_name"]=="eat-2"]["name_1006"].tolist()

bmdf.loc[bmdf["external_gene_name"]=="eat-2", "name_1006"].values[0]

bmdf.to_csv("celegans.go.tsv", sep="\t", index=False)
bmdf.to_excel("celegans.go.xlsx", index=False) # !!! Microsoft Excel has a character limit of 32,767 characters in each cell.

# # DAVID
#
# The Database for Annotation, Visualization and Integrated Discovery (DAVID) provides a comprehensive set of functional annotation tools for investigators to understand the biological meaning behind large lists of genes. These tools are powered by the comprehensive DAVID Knowledgebase built upon the DAVID Gene concept which pulls together multiple sources of functional annotations. For any given gene list, DAVID tools are able to:
#
# - Identify enriched biological themes, particularly GO terms
# - Discover enriched functional-related gene groups
# - Cluster redundant annotation terms
# - Visualize genes on BioCarta & KEGG pathway maps
# - Display related many-genes-to-many-terms on 2-D view.
# - Search for other functionally related genes not in the list
# - List interacting proteins
# - Explore gene names in batch
# - Link gene-disease associations
# - Highlight protein functional domains and motifs
# - Redirect to related literatures
# - Convert gene identifiers from one type to another.
# - And more
#
# [https://david.ncifcrf.gov](https://david.ncifcrf.gov)
#
# [DAVID Web Service](https://david.ncifcrf.gov/content.jsp?file=WS.html)
#
# [DAVID registration](https://david.ncifcrf.gov/webservice/services/DAVIDWebService)
#
# [AGEpy DAVIDenrich function](https://github.com/mpg-age-bioinformatics/AGEpy/blob/663dcce6aecc210a5b57c7e7a3bf66be807eb165/AGEpy/david.py#L19)

import AGEpy as age
import pandas as pd
import random

bmdf=pd.read_csv("celegans.go.tsv", sep="\t")
bmdf.head()

pcgenes=bmdf[bmdf["gene_biotype"]=="protein_coding"]["ensembl_gene_id"].tolist()
query_genes=random.choices(pcgenes, k=250)

help(age.DAVIDenrich)

# {{{
categories = [ 'GOTERM_BP_FAT', 'GOTERM_CC_FAT', 'GOTERM_MF_FAT', \
              'KEGG_PATHWAY', 'BIOCARTA', 'PFAM', 'PROSITE' ]
categories = ",".join(categories)

registered_email="jorge.boucas@age.mpg.de"

david=age.DAVIDenrich('WORMBASE_GENE_ID', categories, registered_email, query_genes, verbose=True)
# }}}

david.head()

# {{{
## add gene names to david output

ids2names=bmdf[["ensembl_gene_id", "external_gene_name"]].astype(str)
ids2names["ensembl_gene_id"]=ids2names["ensembl_gene_id"].apply( lambda x: x.upper() )
ids2names.index=ids2names["ensembl_gene_id"].tolist()
ids2names=ids2names[["external_gene_name"]].to_dict()["external_gene_name"]

def get_names(x,ids2names=ids2names):
    x=x.split(", ")
    x=[ ids2names[s] for s in x ]
    x=", ".join(x)
    return x

david["gene names"]=david["geneIds"].apply( lambda x: get_names(x) )
david.head()
# }}}

# # Image analysis

# {{{
from imutils import contours
from skimage import measure
import numpy as np
import argparse
import imutils
import cv2
import os
import matplotlib.pyplot as plt
import random
import pandas as pd

from scipy.interpolate import splrep, BSpline
from patsy import cr
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# %matplotlib inline
# }}}

# download tif file
import urllib.request 
url="https://datashare.mpcdf.mpg.de/s/mGvEul0YQ4KTtON/download"
IMAGE_FILE = "Cell1.tif"
urllib.request.urlretrieve(url, IMAGE_FILE)

_, images=cv2.imreadmulti(IMAGE_FILE,[],cv2.IMREAD_ANYDEPTH | cv2.IMREAD_UNCHANGED)
print(len(images))

# {{{
# get an image from the series as a test image
image=images[0]

image_original=image.copy()
image_copy = image.copy()

plt.imshow(image,cmap = 'gray', interpolation='bicubic')
plt.title(IMAGE_FILE)
plt.show()
# }}}

# {{{
# # note required in this example
# # convert to gray scale
# image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
# plt.imshow(image,cmap = 'gray', interpolation='bicubic')
# plt.title("gray")
# plt.show()
# }}}

# blur image
image = cv2.GaussianBlur(image, (5, 5), 0)
blur = image.copy()
plt.imshow(image, cmap = 'gray', interpolation='bicubic')
plt.title("blurred")
plt.show()

# threshold
threshold=30
ret3,image = cv2.threshold(image, threshold, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
plt.imshow(image,cmap = 'gray', interpolation='bicubic')
plt.title("thresh")
plt.show()

# erode
erode=2
image = cv2.erode(image, None, iterations=erode)
plt.imshow(image,cmap = 'gray', interpolation='bicubic')
plt.title("erode")
plt.show()

#dilate
dilate=10
image = cv2.dilate(image, None, iterations=dilate)
plt.imshow(image,cmap = 'gray', interpolation='bicubic')
plt.title("dilate")
plt.show()

# {{{
# find contours
image=image.astype(np.uint8)
countours, hierarchy = cv2.findContours(image=image, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_SIMPLE)

## find the biggest counter in image
area= max([ cv2.contourArea(cnt) for cnt in countours ])
biggest_countour = [ cnt for cnt in countours if cv2.contourArea(cnt) == area ]

print(area)

cv2.drawContours(image=image_copy, contours=biggest_countour, contourIdx=-1, color=(0, 255, 0), thickness=1, lineType=cv2.LINE_AA)

plt.imshow(image_copy,cmap = 'gray', interpolation='bicubic')
plt.title(f"{IMAGE_FILE}, area={area}")
plt.show()
# }}}

# {{{
# first frames should be discarded, relevant frames start at 4

plt.imshow(images[3],cmap = 'gray', interpolation='bicubic')
plt.title("first relevant frame")
plt.show()
# }}}

# {{{
df=[]
_, images=cv2.imreadmulti(IMAGE_FILE,[],cv2.IMREAD_ANYDEPTH | cv2.IMREAD_UNCHANGED)
# print(images)
i=3
pos=1
for image in images[i:]:
    i=i+1
    image_original=image.copy()
    image_copy = image.copy()

    # blur image
    image = cv2.GaussianBlur(image, (5, 5), 0)
    blur = image.copy()

    # threshold
    ret3,image = cv2.threshold(image, threshold, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)

    # erode
    image = cv2.erode(image, None, iterations=erode)

    #dilate
    image = cv2.dilate(image, None, iterations=dilate)


    # find contours
    image=image.astype(np.uint8)
    countours, hierarchy = cv2.findContours(image=image, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_SIMPLE)
    # contours, hierarchy = cv2.findContours(image=image, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE)

    ## find the biggest counter in image
    area= max([ cv2.contourArea(cnt) for cnt in countours ])
    biggest_countour = [ cnt for cnt in countours if cv2.contourArea(cnt) == area ]

    df.append([i,area])

    cv2.drawContours(image=image_copy, contours=biggest_countour, contourIdx=-1, color=(0, 255, 0), thickness=1, lineType=cv2.LINE_AA)

df=pd.DataFrame(df,columns=["frame","area"])

df.head()
# }}}

# {{{
m=df["area"].mean()
s=df["area"].std()

df["z-score"]=df["area"].apply(lambda x: (x-m)/s )

x=df["frame"].tolist()
y=df["area"].tolist()
y_=df["z-score"].tolist()

def plot_smoothed(x, y, dof=5, label=None):

    # Generate spline basis with different degrees of freedom
    x_basis = cr(x, df=dof, constraints="center")

    # Fit model to the data
    model = LinearRegression().fit(x_basis, y)

    # Get estimates
    y_hat = model.predict(x_basis)

    plt.plot(x, y_hat, label=label)

i=1
figsize=[12,3]
fig=plt.figure(figsize=(figsize[0], figsize[1]))

plt.subplot(1,2,i)
plt.scatter(df["frame"],df["area"], color="k", s=1)
plot_smoothed(x,y, dof=13)
plt.ylabel("area")
plt.xlabel("frame")
plt.title(f"{IMAGE_FILE} (area)")

i=i+1
plt.subplot(1,2,i)
plt.scatter(df["frame"],df["z-score"], color="k", s=1)
plot_smoothed(x,y_, dof=13)
plt.hlines(0, min(x), max(x), colors="k", linestyles='dashed')
plt.ylabel("z-score")
plt.xlabel("frame")
plt.title(f"{IMAGE_FILE} (z-score)")

plt.tight_layout()
plt.show()
# }}}
# ## hands-on large data: exploring GTEX RNAseq data, finding gene-gene correlations and age related changes
#
# https://gtexportal.org/
#
# Gene TPM file - https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz - saved under: 
# ```
# /mnt/training-volume/common/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct
# ```


import pandas as pd
import dask.dataframe as dd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# gct="/mnt/training-volume/common/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"

file_stats = os.stat(gct)
print(f'File Size in GigaBytes is {file_stats.st_size / (1024 * 1024 * 1024)}')

df=pd.read_csv(gct, skiprows=2, nrows=10, sep="\t")
df.head()

# {{{
# get all genes
genes=pd.read_csv(gct, skiprows=2, usecols=["Name","Description"], sep="\t")
print("Genes:\n", genes.head())

# get all samples
samples=df.columns.tolist()
samples=[ s for s in samples if s not in ["Name", "Description"] ]
print("\n\nSamples:\n", samples[:10])
# }}}

# {{{
fileheader_size=2
table_heder=1
gene_id="ENSG00000227232.5" # WASH7P
gene_index=genes[genes["Name"]==gene_id].index.tolist()[0]

skiprows=fileheader_size+table_heder+gene_index

tmp=pd.read_csv(gct, skiprows=skiprows, nrows=1, usecols=[0,1], sep="\t", header=None)
tmp.head()
# }}}

tmp=pd.read_csv(gct, skiprows=skiprows, nrows=1, names=["Name", "Description"]+samples, sep="\t", header=None)
tmp


def getgene(geneid,samples=samples,genes=genes):
    fileheader_size=2
    table_heder=1
    gene_index=genes[genes["Name"]==geneid].index.tolist()[0]

    skiprows=fileheader_size+table_heder+gene_index
    tmp=pd.read_csv(gct, skiprows=skiprows, nrows=1, names=["Name", "Description"]+samples, sep="\t", header=None)
    return tmp
WASH7P=getgene("ENSG00000227232.5")
WASH7P

sns.kdeplot(WASH7P[samples].transpose()[0].tolist())
plt.title("linear")
plt.show()
sns.kdeplot(np.log10( WASH7P[samples].transpose()[0].tolist()) )
plt.title("log10")
plt.show()


# ## Pearson correlation coefficient
#
# In statistics, the Pearson correlation coefficient (PCC) is a correlation coefficient that measures linear correlation between two sets of data.
#
# ## Spearman's rank correlation coefficient
#
# In statistics, Spearman's rank correlation coefficient or Spearman's ρ, is a nonparametric measure of rank correlation (statistical dependence between the rankings of two variables). 

from scipy.stats import pearsonr 
from scipy.stats import spearmanr
import pandas as pd
import numpy as np

WASH7P=getgene("ENSG00000227232.5")
FAM138A=getgene("ENSG00000237613.2")
pair=pd.concat([WASH7P,FAM138A])
pair.reset_index(inplace=True, drop=True)
pair=pair[samples].transpose()
pair=pair[ ( pair[0] > 0 ) & ( pair[1] > 0 ) ] 
pair.head()

# {{{
print("linear")
print(pearsonr( pair[0].tolist(), pair[1].tolist() ))
print(spearmanr( pair[0].tolist(), pair[1].tolist() ))

print("log10")
pair=np.log2(pair)
print(pearsonr( pair[0].tolist(), pair[1].tolist() ))
print(spearmanr( pair[0].tolist(), pair[1].tolist() ))
# }}}

# {{{
gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct_"
donepairs=[]
i=0
results=[]
with open(gct, "r") as f:
    with open(gct, "r") as f_ :
        for ll in f:
            i=i+1
            i_=0
            if i >= 4:
                for ll_ in f_:
                    i_=i_+1
                    if i_ >= 4:
                        if ( i == i_ ) or ( f"{i}:{i_}" in donepairs ) or ( f"{i_}:{i}" in donepairs ) :
                            continue
                        i_=i_+1
                        
                        l=ll[:-2].split("\t")
                        l_=ll_[:-2].split("\t")
                        gid=l[0]
                        gid_=l_[0]
                        l=l[2:]
                        l_=l_[2:]
                        
                        tmp=pd.DataFrame( {0:l,1:l_} )
                        tmp=tmp=tmp[ ( tmp[0]!="" ) & (tmp[1]!="") ]
                        #print(tmp)
                        tmp=tmp.astype(float)
                        tmp=tmp[ ( tmp[0]>0 ) & (tmp[1]>0) ]
                        
                        n=len(tmp)
                        
                        if n > 2 :
                        
                            tmp=np.log10(tmp)

                            l=tmp[0].tolist()
                            l_=tmp[1].tolist()

                            pearson_stat, pearson_p=pearsonr( l, l_ )
                            spearman_corr, spearman_p=spearmanr( l, l_ )

                            results.append( [ gid, gid_, n, pearson_stat, pearson_p, spearman_corr, spearman_p] )
                        
                        donepairs=donepairs+[ f"{i}:{i_}", f"{i_}:{i}" ]
                        
                        # we are breaking the loop if we are over the 10 line just for demo purposes
                        if i_ > 10:
                        
                            break
                              
results=pd.DataFrame(results, columns=["gene 1", "gene 2", "n", "pearson_stat", "pearson_p", "spearman_corr", "spearman_p"] )  
results
# }}}
# {{{
gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct_"
donepairs=[]
i=0
results=[]
with open(gct, "r") as f:
    with open(gct, "r") as f_ :
        for ll in f:
            i=i+1
            i_=0
            if i >= 4:
                for ll_ in f_:
                    i_=i_+1
                    if i_ >= 4:
                        if ( i == i_ ) or ( f"{i}:{i_}" in donepairs ) or ( f"{i_}:{i}" in donepairs ) :
                            continue
                        i_=i_+1
                        
                        l=ll[:-2].split("\t")
                        l_=ll_[:-2].split("\t")
                        gid=l[0]
                        gid_=l_[0]
                        l=l[2:]
                        l_=l_[2:]
                        
                        tmp=pd.DataFrame( {0:l,1:l_} )
                        tmp=tmp=tmp[ ( tmp[0]!="" ) & (tmp[1]!="") ]
                        #print(tmp)
                        tmp=tmp.astype(float)
                        tmp=tmp[ ( tmp[0]>0 ) & (tmp[1]>0) ]
                        
                        n=len(tmp)
                        
                        if n > 2 :
                        
                            tmp=np.log10(tmp)

                            l=tmp[0].tolist()
                            l_=tmp[1].tolist()

                            pearson_stat, pearson_p=pearsonr( l, l_ )
                            spearman_corr, spearman_p=spearmanr( l, l_ )

                            results.append( [ gid, gid_, n, pearson_stat, pearson_p, spearman_corr, spearman_p] )
                        
                        donepairs=donepairs+[ f"{i}:{i_}", f"{i_}:{i}" ]
                        
                        # we are breaking the loop if we are over the 10 line just for demo purposes
                        if i_ > 10:
                        
                            break
                              
results=pd.DataFrame(results, columns=["gene 1", "gene 2", "n", "pearson_stat", "pearson_p", "spearman_corr", "spearman_p"] )  
results
# }}}

# {{{
import multiprocessing as mp


def blocks(f, cut, size=64*1024): # 65536
    start, chunk =cut
    iter=0
    read_size=int(size)
    _break =False
    while not _break:
        if _break: break
        if f.tell()+size>start+chunk:
            read_size=int(start+chunk- f.tell() )
            _break=True
        b = f.read(read_size)
        iter +=1
        if not b: break
        yield b


def get_chunk_line_count(data):
    fn,  chunk_id, cut = data
    start, chunk =cut
    cnt =0
    last_bl=None

    with open(fn, "r") as f:
        if 0:
            f.seek(start)
            bl = f.read(chunk)
            cnt= bl.count('\n')
        else:
            f.seek(start)
            for i, bl  in enumerate(blocks(f,cut)):
                cnt +=  bl.count('\n')
                last_bl=bl

        if not last_bl.endswith('\n'):
            cnt -=1

        return cnt
....
pool = multiprocessing.Pool(processes=pool_size,
                            initializer=start_process,
                            )
pool_outputs = pool.map(get_chunk_line_count, inputs)
pool.close() # no more tasks
pool.join() 
# }}}

# {{{ jupyter={"outputs_hidden": true}
"""randsamp - extract a random subset of n lines from a large file"""

import random

gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct_"

def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with open(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def sample_lines(path, linepos, nsamp):
    """return nsamp lines from path where line offsets are in linepos"""
    offsets = random.sample(linepos, nsamp)
    print(offsets)
    offsets.sort()  # this may make file reads more efficient

    lines = []
    with open(path) as inf:
        for offset in offsets:
            inf.seek(offset)
            lines.append(inf.readline())
    return lines

dataset = 'big_data.txt'
nsamp = 5
linepos = scan_linepos(gct) # the scan only need be done once

lines = sample_lines(gct, linepos, nsamp)
print('selecting %d lines from a file of %d' % (nsamp, len(linepos)) )
print(''.join(lines) )
# }}}

# {{{
import itertools

gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct_"

n_processors=4

def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with open(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def sample_lines( linepos, path=gct):
    """return nsamp lines from path where line offsets are in linepos"""
    # offsets = random.sample(linepos, nsamp)
    # print(offsets)
    linepos.sort()  # this may make file reads more efficient

    results = []
    offset0=None
    offset1=None
    
    with open(path) as inf:
        
        for offset in linepos:
            
            if offset[0] != offset0 :
                offset0=offset[0]
                inf.seek(offset0)
                l0=inf.readline().split("\n")[0]
                
            if offset[1] != offset1 :
                offset1=offset[1]
                inf.seek(offset1)
                l1=inf.readline().split("\n")[0]
            
            l=l0.split("\t")
            l_=l1.split("\t")
            gid=l[0]
            gid_=l_[0]
            l=l[2:]
            l_=l_[2:]

            tmp=pd.DataFrame( {0:l,1:l_} )
            tmp=tmp[ ( tmp[0]!="" ) & (tmp[1]!="") ]
            #print(tmp)
            tmp=tmp.astype(float)
            tmp=tmp[ ( tmp[0]>0 ) & (tmp[1]>0) ]

            n=len(tmp)

            if n > 2 :

                tmp=np.log10(tmp)

                l=tmp[0].tolist()
                l_=tmp[1].tolist()

                pearson_stat, pearson_p=pearsonr( l, l_ )
                spearman_corr, spearman_p=spearmanr( l, l_ )
                
                res=[ gid, gid_, n, pearson_stat, pearson_p, spearman_corr, spearman_p]
                res=[ str(s) for s in res ] 
                res="\t".join(res)

                results.append( res )
            
            # lines.append(inf.readline())
    return "\n".join(results)

linepos = scan_linepos(gct) # the scan only need be done once
linepos = linepos[3:]
combinations = list(itertools.combinations(linepos, 2))
print(len(combinations))

# combinations=combinations[:40]

lines_to_process=[combinations[i:i + n_processors] for i in range(0, len(combinations), n_processors)]

# print(lines_to_process)

# results = sample_lines(lines_to_process)


pool = mp.Pool(n_processors)
results = []
for d in lines_to_process:
    output = pool.apply_async(sample_lines, [d])
    results.append(output)
pool.close() # no more tasks
# pool.join()

results=[ s.get() for s in results ] 
results="\n".join(results)

print(results)
# }}}

pd.DataFrame([(1,2),(3,4)])

# {{{
import itertools

gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct_"

n_processors=4

def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with open(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def sample_lines( linepos, path=gct):
    """return nsamp lines from path where line offsets are in linepos"""

    linepos.sort()  # this may make file reads more efficient

    results = []
    
    linepos=pd.DataFrame(linepos)
    
    for offset0 in list(set(linepos[0].tolist() )):
        
        with open(path) as inf:
            inf.seek(offset0)
            l0=inf.readline().split("\n")[0]
        
        def _corr(offset1,l0=l0):
            
            with open(path) as inf:

                inf.seek(offset1)
                l1=inf.readline().split("\n")[0]

                l=l0.split("\t")
                l_=l1.split("\t")
                gid=l[0]
                gid_=l_[0]
                l=l[2:]
                l_=l_[2:]

                tmp=pd.DataFrame( {0:l,1:l_} )
                tmp=tmp[ ( tmp[0]!="" ) & (tmp[1]!="") ]
                #print(tmp)
                tmp=tmp.astype(float)
                tmp=tmp[ ( tmp[0]>0 ) & (tmp[1]>0) ]

                n=len(tmp)

                if n > 2 :

                    tmp=np.log10(tmp)

                    l=tmp[0].tolist()
                    l_=tmp[1].tolist()

                    pearson_stat, pearson_p=pearsonr( l, l_ )
                    spearman_corr, spearman_p=spearmanr( l, l_ )

                    res=[ gid, gid_, n, pearson_stat, pearson_p, spearman_corr, spearman_p]


                else:
                    res=[ gid, gid_, n, None, None, None, None]

                res=[ str(s) for s in res ] 
                res="\t".join(res)
                
            return res
        
        linepos_=linepos[linepos[0]==offset0]
        linepos_[2]=linepos[1].apply(lambda x: _corr(x) )
        
        r0="\n".join( linepos_[2].tolist() )
        
        results.append(r0)     
                        
    return "\n".join(results)

linepos = scan_linepos(gct) # the scan only need be done once
linepos = linepos[3:] # remove the header lines
combinations = list(itertools.combinations(linepos, 2))
print(len(combinations), "combinations")

# dev/demo only do the first 40 combinations
combinations=combinations[:40]

# we want to have chuncks of len 500 each
target=50

lines_to_process=[combinations[i:i + target] for i in range(0, len(combinations), target)]

lines_to_process=[lines_to_process[i:i + n_processors] for i in range(0, len(lines_to_process), n_processors)]


results = []

for lines_to_process_ in lines_to_process :
                
    pool = mp.Pool(n_processors)
    for d in lines_to_process_:
        output = pool.apply_async(sample_lines, [d])
        # results_=[ s.get() for s in output ] 
        results.append(output)
    pool.close() # no more tasks
    # pool.join()

results=[ s.get() for s in results ] 
results="\n".join(results)

results=results.split("\n")
results=[ s.split("\t") for s in results ]
results=pd.DataFrame(results, columns=[ "gid", "gid_", "n", "pearson_stat", "pearson_p", "spearman_corr", "spearman_p"])
results.to_csv("gtex.corr.tsv", index=None, sep="\t")
results.to_excel("gtex.corr.xlsx", index=None)
results.head()
# }}}

import os
os.path.isfile("/")

# {{{
# %%writefile gtex.corr.py
import warnings
warnings.filterwarnings("ignore")
from scipy.stats import pearsonr 
from scipy.stats import spearmanr
from datetime import datetime
from time import process_time 
import multiprocessing as mp
import pandas as pd
import numpy as np
import itertools
import sys
import os
import time

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string, ":: started")

# Start the stopwatch / counter  
t1_start = process_time()  

gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
output="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/gtexcorr/"
if not os.path.isdir( output ):
    os.makedirs( output )

n_processors=32

def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with open(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def sample_lines( linepos_input, p, path=gct):
    """return nsamp lines from path where line offsets are in linepos"""
    
    target=5000

    linepos_input=[linepos_input[i:i + target] for i in range(0, len(linepos_input), target)]
    
    # resulst_=[]
    c=0
    for linepos in linepos_input:
        c=c+target

        linepos.sort()  # this may make file reads more efficient

        filename=f"{output}{linepos[0][0]}.{linepos[0][1]}.tsv"

        if os.path.isfile(filename) :
            return ""

        results = []

        linepos=pd.DataFrame(linepos)

        for offset0 in list(set(linepos[0].tolist() )):

            with open(path) as inf:
                inf.seek(offset0)
                l0=inf.readline().split("\n")[0]

            def _corr(offset1,l0=l0):

                with open(path) as inf:

                    inf.seek(offset1)
                    l1=inf.readline().split("\n")[0]

                    l=l0.split("\t")
                    l_=l1.split("\t")
                    gid=l[0]
                    gid_=l_[0]
                    l=l[2:]
                    l_=l_[2:]

                    tmp=pd.DataFrame( {0:l,1:l_} )
                    tmp=tmp[ ( tmp[0]!="" ) & (tmp[1]!="") ]
                    #print(tmp)
                    tmp=tmp.astype(float)
                    tmp=tmp[ ( tmp[0]>0 ) & (tmp[1]>0) ]

                    n=len(tmp)

                    if n > 2 :

                        tmp=np.log10(tmp)

                        l=tmp[0].tolist()
                        l_=tmp[1].tolist()

                        pearson_stat, pearson_p=pearsonr( l, l_ )
                        spearman_corr, spearman_p=spearmanr( l, l_ )

                        res=[ gid, gid_, n, pearson_stat, pearson_p, spearman_corr, spearman_p]

                    else:
                        res=[ gid, gid_, n, None, None, None, None]

                    res=[ str(s) for s in res ] 
                    res="\t".join(res)

                return res

            linepos_=linepos[linepos[0]==offset0]
            linepos_[2]=linepos_[1].apply(lambda x: _corr(x) )

            r0="\n".join( linepos_[2].tolist() )

            results.append(r0)
        
        results="\n".join(results)
        
        # results_.append(results)

        with open(filename, "w") as f:
            f.write(results)
        
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", "part:", p,",", c, "done" )
        sys.stdout.flush()
        
        
    # results_="\n".join(results_)
                        
    return p

linepos = scan_linepos(gct) # the scan only need be done once
linepos = linepos[3:] # remove the header lines
combinations = list(itertools.combinations(linepos, 2))
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", len(combinations), "combinations" )
sys.stdout.flush()

# dev/demo only do the first 40 combinations
# combinations=combinations[:40]

# we want to have chuncks of len 500 each
# target=5000

# lines_to_process=[combinations[i:i + target] for i in range(0, len(combinations), target)]

size_blocks=int(len(combinations)/n_processors)

lines_to_process=[combinations[i:i + size_blocks] for i in range(0, len(combinations), size_blocks)]

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", len(lines_to_process), "parts" )
sys.stdout.flush()

results = []

i=0      

# for lines_to_process_ in lines_to_process :
    
    # i=i + len(lines_to_process_)*target

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", "starting pool" )
sys.stdout.flush()
    
pool = mp.Pool(n_processors)
for d in lines_to_process:
    output = pool.apply_async(sample_lines, [d,i])
    i=i+1
    results.append(output)
pool.close() # no more tasks
time.sleep(10)
    # print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", i, "done" )
    # sys.stdout.flush()
    # pool.join()

results=[ s.get() for s in results ] 
results="\n".join(results)
# print(f"Done\n:{results}"

# results=results.split("\n")
# results=[ s.split("\t") for s in results ]
# results=pd.DataFrame(results, columns=[ "gid", "gid_", "n", "pearson_stat", "pearson_p", "spearman_corr", "spearman_p"])
# results.to_csv("gtex.corr.all.tsv", index=None, sep="\t")
# results.to_excel("gtex.corr.all.xlsx", index=None)

# Stop the stopwatch / counter 
t1_stop = process_time() 

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string, ":: finished", )
print("Elapsed time:", t1_stop, t1_start) 
sys.exit(0)
# }}}

# %run gtex.corr.py

df = dd.read_csv(gct, sample=2560000,sep="\t", skiprows=2)

df.head()

combinations = list(itertools.combinations(linepos, 2))
combinations[0][0]

pool.close()


