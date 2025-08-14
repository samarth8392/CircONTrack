# Coverage Peak Detection - Explained for Non-Statisticians

## What This Code Does
This tool finds "peaks" in DNA sequencing data - regions where there's unusually high read coverage that might indicate circular DNA, duplications, or other interesting genomic features.

## Key Libraries and Their Roles

- **pysam**: Reads BAM files (compressed sequencing alignment files)
- **scipy.stats**: Provides statistical distributions and tests
- **numpy**: Mathematical operations and array handling
- **statsmodels**: Multiple testing correction (prevents false discoveries)
- **matplotlib/seaborn**: Optional plotting capabilities

---

## Step 1: Understanding the Data Distribution Problem

### The Basic Concept
```python
# If data variance > mean → Negative Binomial distribution (handles noisy/overdispersed data)
def fit_negative_binomial(self, coverage_values):
    coverage_values = np.array(coverage_values)
    non_zero = coverage_values[coverage_values > 0]
    
    mean_cov = np.mean(non_zero)      # Average coverage
    var_cov = np.var(non_zero)        # How spread out the data is
    
    # The key check: Is variance bigger than mean?
    if var_cov <= mean_cov:
        # Data is "well-behaved" - use simpler Poisson distribution
        return {'distribution': 'poisson', 'mean': mean_cov}
    
    # Data is "overdispersed" - use Negative Binomial
    p = mean_cov / var_cov
    n = mean_cov * p / (1 - p)
```

**What "overdispersed" means:**
- **Normal data**: If you flip 100 coins, you expect ~50 heads with small variation
- **Overdispersed data**: Like sequencing - some regions get way more reads than expected by chance
- **Why it happens**: PCR bias, repetitive sequences, technical artifacts

---

## Step 2: Sliding Window Coverage Calculation

### How Coverage is Calculated
```python
def calculate_coverage_windows(self, chromosome, start=None, end=None):
    """
    Think of this like scanning a chromosome with a magnifying glass:
    - Move a "window" (e.g., 1000bp) across the chromosome
    - Count how many reads overlap each window
    - Calculate average coverage for that window
    """
    windows = []
    
    # Slide the window across the chromosome
    for window_start in range(start, end, self.step_size):
        window_end = min(window_start + self.window_size, end)
        
        read_count = 0
        coverage_sum = 0
        
        # Count reads that overlap this window
        for read in self.bam.fetch(chromosome, window_start, window_end):
            if read.is_unmapped or read.is_secondary:
                continue  # Skip bad reads
                
            # Calculate how much of the read overlaps our window
            read_start = max(read.reference_start, window_start)
            read_end = min(read.reference_end, window_end)
            overlap = read_end - read_start
            
            if overlap > 0:
                read_count += 1
                coverage_sum += overlap
        
        # Average coverage = total overlapping bases / window size
        avg_coverage = coverage_sum / window_length
```

**What this creates:**
- A list of windows, each with its coverage value
- Example: Window 1: 5.2X, Window 2: 4.8X, Window 3: 15.7X ← This might be a peak!

---

## Step 3: Fitting the Statistical Background Model

### What "Fitting" Means
```python
# Method of moments estimation - a simple way to fit the distribution
mean_cov = np.mean(non_zero)     # What's the typical coverage?
var_cov = np.var(non_zero)       # How much does it vary?

# Calculate Negative Binomial parameters
# These formulas come from statistics theory
p = mean_cov / var_cov           # "Success probability" 
n = mean_cov * p / (1 - p)       # "Number of failures"
```

**Think of it like this:**
1. **Background model**: "What does normal coverage look like across the genome?"
2. **Parameters**: Numbers that describe the typical pattern
3. **Fitting**: Finding the best numbers that match your data

**Real example:**
- Your data: Most windows have 8-12X coverage, but some have 2X, others have 40X
- The model learns: "Typical coverage is ~10X, but there's lots of variation"
- Uses this to decide: "40X coverage is way higher than expected!"

---

## Step 4: Statistical Significance Testing

### P-value Calculation
```python
def calculate_pvalues(self, windows, dist_params):
    """
    For each window, ask: "What's the probability of seeing this much
    coverage (or higher) if this was just normal background?"
    """
    for window in windows:
        coverage = window['coverage']
        
        if dist_params['distribution'] == 'nbinom':
            # "What's the chance of getting ≥ this coverage by random chance?"
            pval = stats.nbinom.sf(coverage - 1, 
                                  dist_params['n'], 
                                  dist_params['p'])
```

**P-value interpretation:**
- **p = 0.001**: Only 0.1% chance this is random → Probably a real peak!
- **p = 0.3**: 30% chance this is random → Probably just noise
- **p = 0.05**: The typical cutoff (5% chance of being wrong)

---

## Step 5: Multiple Testing Correction

### Why We Need This
```python
# Problem: If you test 10,000 windows at p < 0.05, you'll get ~500 false positives!
# Solution: Adjust p-values to control false discovery rate

rejected, adjusted_pvals, _, _ = multipletests(
    pvalues, 
    alpha=fdr_threshold,  # Usually 0.05 = 5% false discovery rate
    method='fdr_bh'       # Benjamini-Hochberg correction
)
```

**The multiple testing problem:**
- Test 1 window: 5% chance of false positive
- Test 10,000 windows: You'll get hundreds of false positives!
- **Solution**: Make the cutoff stricter when testing many things

---

## Step 6: Peak Filtering and Merging

### Applying Biological Filters
```python
# Don't just rely on statistics - apply biological sense
for window in windows:
    fold_change = window['coverage'] / background
    
    # Keep peaks that are:
    if (is_statistically_significant and          # p < 0.05 (adjusted)
        fold_change >= min_fold_change and        # At least 2x higher than background  
        window['coverage'] >= min_coverage):      # At least 5X coverage (not just noise)
        
        peaks.append(window)
```

### Merging Nearby Peaks
```python
# If two peaks are close together (< 1000bp), merge them
# This prevents splitting one real peak into multiple small ones

if peak['start'] - current['end'] <= max_distance:
    # Merge: extend the peak boundaries
    current['end'] = peak['end']
    current['coverage'] = max(current['coverage'], peak['coverage'])
```

---

## Summary: The Complete Pipeline

1. **Slide windows** across chromosomes, calculate coverage
2. **Fit statistical model** to learn what "normal" coverage looks like
3. **Calculate p-values** for each window: "Is this coverage unusual?"
4. **Correct for multiple testing** to avoid too many false positives
5. **Apply filters**: Must be statistically significant AND biologically meaningful
6. **Merge nearby peaks** to get final peak calls
7. **Output BED file** with peak locations and statistics

## Key Insight
This isn't just finding "high coverage" - it's finding "**statistically unexpected** high coverage" while being careful about false discoveries. The negative binomial distribution helps handle the messy, noisy nature of real sequencing data.