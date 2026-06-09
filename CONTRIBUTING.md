# Contributing to ClarVar

Thank you for your interest in contributing to ClarVar! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)
- [Types of Contributions](#types-of-contributions)

---

## Code of Conduct

We are committed to providing a welcoming and inspiring community for all. Please be respectful and constructive in all interactions.

### Our Pledge

- We welcome contributions from people of all backgrounds and experience levels
- We are committed to creating a safe, inclusive community
- We will address concerns fairly and promptly

---

## Getting Started

### Prerequisites

- Python 3.8+
- Git
- GitHub account
- Basic knowledge of Python

### Development Setup

```bash
# 1. Fork the repository on GitHub
#    (Click "Fork" button on https://github.com/clarvar/clarvar)

# 2. Clone your fork locally
git clone https://github.com/YOUR-USERNAME/clarvar.git
cd clarvar

# 3. Create a virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 4. Install in development mode with dependencies
pip install -e ".[dev]"

# 5. Verify setup
python -m pytest tests/
```

---

## Development Workflow

### 1. Create a Feature Branch

```bash
# Make sure you're on main and up to date
git checkout main
git pull origin main

# Create a feature branch
git checkout -b feature/my-feature
```

**Branch naming conventions:**
- Features: `feature/description`
- Bug fixes: `bugfix/description`
- Documentation: `docs/description`
- Tests: `test/description`

### 2. Make Your Changes

Write code, tests, and documentation as needed.

### 3. Run Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=clarvar tests/

# Run specific test file
pytest tests/test_variant.py

# Run with verbose output
pytest -v
```

### 4. Check Code Style

```bash
# Format code with black
black src/clarvar tests/

# Check with flake8
flake8 src/clarvar tests/

# Type checking with mypy
mypy src/clarvar
```

### 5. Commit Your Changes

```bash
git add .
git commit -m "Add feature: description"
```

**Commit message guidelines:**
- Start with a verb: "Add", "Fix", "Update", "Refactor"
- Be descriptive but concise
- Reference issues: "Fixes #123"
- Examples:
  - "Add gnomAD API integration"
  - "Fix variant parsing for multi-allelic sites"
  - "Update documentation for prioritization"

### 6. Push and Create Pull Request

```bash
# Push to your fork
git push origin feature/my-feature

# Then create a Pull Request on GitHub
```

---

## Coding Standards

### Python Style

We follow PEP 8 with these guidelines:

```python
# Good
def calculate_priority_score(variant: Variant, weights: Dict[str, float]) -> float:
    """
    Calculate priority score for a variant.
    
    Args:
        variant: Variant to score
        weights: Dictionary of weight parameters
        
    Returns:
        Priority score (0-100)
    """
    consequence_score = calculate_consequence_score(variant)
    return consequence_score * weights.get('consequence', 0.4)


# Bad
def calc_score(v, w):  # unclear naming
    return v.cs * w['c']  # cryptic abbreviations
```

### Documentation Strings

Use Google-style docstrings:

```python
def my_function(param1: str, param2: int) -> bool:
    """
    Brief description of what the function does.
    
    Longer description if needed. Explain the purpose,
    behavior, and any important details.
    
    Args:
        param1: Description of param1
        param2: Description of param2
        
    Returns:
        Description of return value
        
    Raises:
        ValueError: Description of when this is raised
        
    Example:
        >>> result = my_function("test", 42)
        >>> print(result)
        True
    """
    pass
```

### Type Hints

Always use type hints:

```python
from typing import List, Optional, Dict

def process_variants(variants: List[Variant]) -> Dict[str, int]:
    """Process variants and return summary stats."""
    pass

def get_gene(variant: Optional[Variant]) -> Optional[str]:
    """Get gene name from variant."""
    pass
```

### Naming Conventions

```python
# Variables and functions: snake_case
variant_list = []
def calculate_score():
    pass

# Classes: PascalCase
class VariantCollection:
    pass

# Constants: UPPER_SNAKE_CASE
DEFAULT_THRESHOLD = 0.01
MAX_VARIANTS = 1000000
```

---

## Testing

### Writing Tests

```python
import pytest
from clarvar.variant import Variant, Consequence

def test_variant_creation():
    """Test creating a variant."""
    variant = Variant(
        chromosome="1",
        position=100,
        ref="A",
        alt="G"
    )
    assert variant.chromosome == "1"

def test_variant_validation():
    """Test variant validation."""
    with pytest.raises(ValueError):
        Variant(chromosome="1", position=0, ref="A", alt="G")

class TestVariantCollection:
    """Test cases for VariantCollection."""
    
    def test_add_variant(self):
        """Test adding a variant to collection."""
        pass
```

### Test Organization

```
tests/
├── test_variant.py       # Tests for variant.py
├── test_annotator.py     # Tests for annotator.py
├── test_prioritizer.py   # Tests for prioritizer.py
├── test_vcf_parser.py    # Tests for vcf_parser.py
├── test_cli.py           # Tests for cli.py
└── fixtures/             # Test data files
    ├── sample.vcf
    └── expected_output.vcf
```

### Test Coverage

Aim for 80%+ code coverage:

```bash
# Generate coverage report
pytest --cov=clarvar --cov-report=html tests/

# View coverage
open htmlcov/index.html
```

---

## Documentation

### Documentation Files

- **README.md**: Project overview and quick start
- **docs/GETTING_STARTED.md**: User guide
- **docs/ARCHITECTURE.md**: Technical design
- **Docstrings**: In-code documentation

### Updating Documentation

1. **For new features**: Add examples and documentation
2. **For changes**: Update relevant docs and examples
3. **API changes**: Update API reference

### Documentation Format

Use Markdown with clear structure:

```markdown
# Main Title

Brief description.

## Section

More details with code examples:

\`\`\`python
from clarvar import Variant
v = Variant(...)
\`\`\`

### Subsection

Additional information.
```

---

## Submitting Changes

### Pull Request Process

1. **Create a PR**: Include a clear description
2. **Title**: Describe what the PR does
3. **Description**: 
   - What does it do?
   - Why is it needed?
   - How does it work?
   - Any breaking changes?

```markdown
## Description

This PR adds gnomAD API integration for real variant annotation.

## Changes

- Implement `GnomadAPI` class
- Add frequency lookup methods
- Cache API responses

## Related Issues

Fixes #123

## Testing

- Added tests for API integration
- Tested with sample VCF file
- Verified error handling

## Checklist

- [x] Tests pass
- [x] Code follows style guidelines
- [x] Documentation updated
- [x] Commits have clear messages
```

### Code Review

All PRs will be reviewed by maintainers. We may request changes:

- **Style issues**: Code format, naming, structure
- **Logic issues**: Correctness, efficiency, edge cases
- **Testing**: Coverage, edge cases, error handling
- **Documentation**: Clarity, completeness

---

## Types of Contributions

### 1. Bug Reports

Found a bug? Create an issue with:
- Description of the problem
- Steps to reproduce
- Expected vs. actual behavior
- Your environment (OS, Python version)

### 2. Feature Requests

Have an idea? Suggest it as an issue:
- Description of the feature
- Use cases and why it's needed
- Proposed implementation (optional)

### 3. Documentation

Improve docs by:
- Fixing typos and unclear sections
- Adding examples
- Creating tutorials
- Improving API documentation

### 4. Code Improvements

- Optimize performance
- Refactor for clarity
- Fix code style issues
- Improve error handling
- Add missing tests

### 5. New Features

Major features should start with an issue for discussion:
- Discuss design and approach
- Get feedback from maintainers
- Plan implementation
- Then submit PR

---

## Project Structure

Understanding the project layout:

```
clarvar/
├── src/clarvar/          # Main package
│   ├── variant.py        # Variant data models
│   ├── annotator.py      # Annotation logic
│   ├── prioritizer.py    # Prioritization algorithms
│   ├── vcf_parser.py     # File parsing
│   ├── cli.py            # Command-line interface
│   └── __init__.py       # Package init
│
├── tests/                # Unit tests
│   ├── test_variant.py
│   ├── test_annotator.py
│   ├── test_prioritizer.py
│   ├── test_vcf_parser.py
│   └── test_cli.py
│
├── docs/                 # Documentation
│   ├── GETTING_STARTED.md
│   └── ARCHITECTURE.md
│
├── examples/             # Example scripts and data
│   ├── sample_variants.vcf
│   └── example_usage.py
│
├── setup.py              # Package configuration
├── requirements.txt      # Dependencies
├── README.md             # Project readme
├── CONTRIBUTING.md       # This file
└── LICENSE               # MIT License
```

---

## Common Tasks for Contributors

### Adding a New Feature

1. Create feature branch: `git checkout -b feature/my-feature`
2. Implement feature with tests
3. Update documentation
4. Commit with clear message: `git commit -m "Add feature: description"`
5. Push and create PR

### Fixing a Bug

1. Create branch: `git checkout -b bugfix/issue-description`
2. Write test that reproduces bug
3. Fix the bug
4. Verify test passes
5. Submit PR

### Improving Documentation

1. Create branch: `git checkout -b docs/improvement`
2. Update relevant docs
3. Verify formatting looks good
4. Submit PR

### Running a Code Review

1. Clone the PR branch
2. Run tests: `pytest`
3. Check code style: `black --check` + `flake8`
4. Review changes manually
5. Provide constructive feedback

---

## Resources

- **Code Style**: https://pep8.org/
- **Google Docstring Style**: https://google.github.io/styleguide/pyguide.html
- **Testing with pytest**: https://docs.pytest.org/
- **Type Hints**: https://docs.python.org/3/library/typing.html
- **Git Workflow**: https://guides.github.com/introduction/flow/

---

## Questions?

- **Documentation**: Check [README.md](README.md) and [docs/](docs/)
- **Issues**: Search existing issues on GitHub
- **Discussions**: Ask on GitHub Discussions
- **Email**: clarvar@example.com

---

## Recognition

We appreciate all contributions! Contributors will be:
- Credited in README
- Mentioned in release notes
- Added to contributors list

---

## Final Notes

- Start small: fix typos, improve docs
- Read existing code to understand style
- Ask questions: no question is too simple
- Be patient: reviews take time
- Have fun: contributing should be enjoyable!

---

**Thank you for contributing to ClarVar! 🧬**

Together, we're building tools that improve clinical genomics.

