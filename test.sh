#!/bin/bash

if test "$1" = "-t"; then
  sage -t monomials.py poly_free_module.py poly_utils.py
else
  sage test_script.py
fi
