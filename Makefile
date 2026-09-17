.PHONY: examples test lint clean

# Run all worked examples
examples:
	python -m examples.automotive_examples

# Run tests (when they exist)
test:
	python -m pytest tests/ -v

# Basic lint check
lint:
	python -m py_compile thermal_mesh_calculators/__init__.py
	python -m py_compile thermal_mesh_calculators/constants.py
	python -m py_compile thermal_mesh_calculators/conduction.py
	python -m py_compile thermal_mesh_calculators/convection.py
	python -m py_compile thermal_mesh_calculators/radiation.py
	python -m py_compile thermal_mesh_calculators/shields.py
	python -m py_compile thermal_mesh_calculators/transient.py
	python -m py_compile examples/automotive_examples.py
	@echo "All modules compile OK"

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
