## In scientific computing, it is very important to know how to write succint math in LaTex. I use Overleaf,

### Official pages for reference:

1. [Mathematical expression (inline and block)](https://www.overleaf.com/learn/latex/Mathematical_expressions)
2. [Greek letter and math symbols (main)](https://www.overleaf.com/learn/latex/List_of_Greek_letters_and_math_symbols)



### Additional notes:

1. paragraph indent and paragraph skip

```
\setlength{\parindent}{0em}
\setlength{\parskip}{1em}
```


2. bold math and upright (not tilted)

```
\mathbf{}
\mathsf{}
```

3. notes

```
\underset{}{}   # underset, the first bracket represent stuff under the second bracket
\lfloor \rfloor   # floor/ceil
\hat     # hat
\Bar     # bar
\left\| \right\|    # norm
\sim    # ~
\arg    # arg
\mathbbm{1}  # indicator function   # make sure to usepackage{bbm}
\cdot   # a dot that is vertically aligned in the center
\backslash  # \
```

4. [space] in math mode(https://www.overleaf.com/learn/latex/Spacing_in_math_mode)



