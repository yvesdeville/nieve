# Changes in version 0.1.5

## Enhancements

- The `pGEV` function now computes the Hessian, and the result is
  tested.

# Changes in version 0.1.4

## Bug fix

- When the argument `x` contained non-finite values `dGEV` and `dexp1`
  returned incorrect values for `x = -Inf` and `x = Inf`.


# Changes in version 0.1.3

## Enhancements

- Non-finite values (hence `NA`s) have been re-allowed in the
  parameters after some checks.
  
  
# Changes in version 0.1.2

## Enhancements

- Non-finite values (hence `NA`s) are no longer allowed in the
  parameters.
	
- The `rGEV` function no longer returns a named vector.

- The probability functions `p`, `d`, `q` propagate the names of their
  first argument to the result and the `"gradient"` and `"hessian"`
  attributes. This may be useful when working with non-stationary
  Extreme-Value models.
  
- The simulation functions `rGEV`, `rGPD2` and `rexp1` have an `array`
  argument to use the vectors of parameters in a fashion that is
  suitable for non-stationary models, then returning a matrix instead
  of a vector.
  
	
## Bug fixes

- The probability `p`, `d`, `q` functions did not work with `deriv =
  TRUE` when their first argument contained `NA`s. This is now
  possible, the corresponding rows of the `"gradient"` attribute and
  the corresponding slices of `"hessian"` attribute contain only NAs.
