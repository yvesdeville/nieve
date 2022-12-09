
# Changes in version 0.1.2

## Enhancements

- Non-finite values (hence `NA`s) are no longer allowed in the
  parameters.
	
- The `rGEV` function no longer returns a named vector.

- The probability functions `p`, `d`, `q` propagate the names of their
  first argument to the result and the `"gradient"` and `"hessian"`
  attributes. This may be useful when working with non-stationary
  Extreme-Value models.
	
## Bug fixes

- The probability `p`, `d`, `q` functions did not work with `deriv =
  TRUE` when their first argument contained `NA`s. This is now
  possible, the corresponding rows of the `"gradient"` attribute and
  the corresponding slices of `"hessian"` attribute contain only NAs.
