.PHONY: examples test doctest lint typecheck clean

# Run all worked examples
examples:
	python -m examples.automotive_examples

# Run the test suite
test:
	python -m pytest tests/ -v

# Run the README's examples as doctests, as CI does
doctest:
	python -m pytest --doctest-glob=README.md README.md

# Lint with ruff, as CI does (the rule set is in pyproject.toml; needs ruff)
lint:
	python -m ruff check .

# Type-check the package with mypy, as CI does (configured in pyproject.toml; needs mypy)
typecheck:
	python -m mypy

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
