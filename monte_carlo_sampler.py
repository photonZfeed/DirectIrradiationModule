from tabnanny import verbose
import numpy as np
from itertools import combinations
import concurrent.futures
import os
import random
import time
from typing import List, Tuple, Optional, Dict, Any
from collections import defaultdict
from utils.grid import Grid
from systematic_sampler import Sampler 
from numba import njit

class MonteCarloSampler():
    """
    Generates random, symmetric LED configurations with variable LED counts and heights using Monte Carlo sampling.

    This class supports reproducible sampling via seeding and enforces point symmetry and quadrant balance on a square grid.
    Configurations are generated such that the spatial distribution of LEDs is symmetric about the grid center, and the number of LEDs per quadrant is balanced.

    The grid is assumed to be square and centered at (0, 0). The number of LEDs and their heights can be varied within user-specified ranges. Sampling can be performed uniformly or weighted by the number of possible configurations for each LED count.

    :param G: Grid object defining the possible LED positions. Must be a square grid with an odd number of points per side.
        (type: Grid)
    :param led_min: Minimum number of LEDs to place in a configuration (inclusive). Must be >= 1. (type: int, default: 1)
    :param led_max: Maximum number of LEDs to place in a configuration (inclusive). Must be >= led_min. (type: int, default: 16)
    :param height_min: Minimum height above the grid for any LED, in centimeters. (type: float, default: 5.0)
    :param height_max: Maximum height above the grid for any LED, in centimeters. Must be >= height_min. (type: float, default: 15.0)
    :param height_step: Step size for possible LED heights, in centimeters. (type: float, default: 1.0)
    :param seed: Random seed for reproducible results. If None, results are non-reproducible. (type: Optional[int], default: None)
    :param sampling_mode: Sampling strategy for LED counts. "uniform" yields roughly equal samples per LED count; "weighted" samples proportional to the number of possible configurations. (type: str, default: "uniform")
    :param verbose: If True, prints initialization and sampling details to stdout. (type: bool, default: True)

    :raises ValueError: If the grid is not square with an odd number of points per side, or if parameter constraints are violated.

    :note:
        - The grid must be square and have an odd number of points per side to ensure a unique center point for symmetry.
        - The class precomputes all valid symmetric configurations for each LED count in the specified range.
        - All random sampling is performed using the provided seed for reproducibility, if specified.
    """
    
    def __init__(self, G: Grid, led_min: int = 1, led_max: int = 16,
                 height_min: float = 5.0, height_max: float = 15.0, 
                 height_step: float = 1, seed: Optional[int] = None,
                 sampling_mode: str = "uniform", verbose: bool = True):
        """
        Constructs a MonteCarloSampler instance for generating symmetric LED configurations on a square grid.

        The sampler enforces point symmetry and quadrant balance, and supports reproducible random sampling. All valid symmetric configurations for each LED count are precomputed at initialization.

        :param G: Grid object defining the possible LED positions. Must be a square grid with an odd number of points per side.
            (type: Grid)
        :param led_min: Minimum number of LEDs to place in a configuration (inclusive). Must be >= 1. (type: int, default: 1)
        :param led_max: Maximum number of LEDs to place in a configuration (inclusive). Must be >= led_min. (type: int, default: 16)
        :param height_min: Minimum height above the grid for any LED, in centimeters. (type: float, default: 5.0)
        :param height_max: Maximum height above the grid for any LED, in centimeters. Must be >= height_min. (type: float, default: 15.0)
        :param height_step: Step size for possible LED heights, in centimeters. (type: float, default: 1.0)
        :param seed: Random seed for reproducible results. If None, results are non-reproducible. (type: Optional[int], default: None)
        :param sampling_mode: Sampling strategy for LED counts. "uniform" yields roughly equal samples per LED count; "weighted" samples proportional to the number of possible configurations. (type: str, default: "uniform")
        :param verbose: If True, prints initialization and sampling details to stdout. (type: bool, default: True)

        :raises ValueError: If the grid is not square with an odd number of points per side, or if parameter constraints are violated.

        :side effects:
            - Prints initialization details to stdout if `verbose` is True.
            - Sets the global NumPy random seed if `seed` is provided.

        :note:
            - The grid must be square and have an odd number of points per side to ensure a unique center point for symmetry.
            - The class precomputes all valid symmetric configurations for each LED count in the specified range.
            - All random sampling is performed using the provided seed for reproducibility, if specified.
        """
        if np.sqrt(len(G.grid_points)).astype(int) % 2 == 0:
            raise ValueError("Grid size must be odd for clear center point")

        if led_min < 1:
            raise ValueError("Minimum LED count must be at least 1")
        if led_max < led_min:
            raise ValueError("Maximum LED count must be greater than or equal to minimum")
        if height_max < height_min:
            raise ValueError("Maximum height must be greater than or equal to minimum height")
        if sampling_mode not in ["uniform", "weighted"]:
            raise ValueError("sampling_mode must be 'uniform' or 'weighted'")
            
        self.G = G
        self.grid_size = np.sqrt(len(G.grid_points)).astype(int)  # Assuming square grid
        self.grid_step = G.step

        self.led_min = led_min
        self.led_max = led_max
        self.height_min = height_min
        self.height_max = height_max
        self.height_step = height_step
        self.seed = seed
        self.sampling_mode = sampling_mode
        self.verbose = verbose
        
        # Initialize random state
        self._random_state = random.Random(seed)
        if seed is not None:
            np.random.seed(seed)
        
        # Generate possible height values
        self.possible_heights = self._generate_height_range()
        self.half_size = self.grid_size // 2
        self.coords = list(range(-self.half_size, self.half_size + 1))
        
        # Precompute all point classifications
        self.Q1_points, self.Q2_points, self.Q3_points, self.Q4_points = [], [], [], []
        self.pos_x_axis, self.pos_y_axis, self.all_pos_axis = [], [], []
        self.neg_x_axis, self.neg_y_axis = [], []
        self.center = (0, 0)  
        self._classify_points()
        
        # Precompute valid combinations for each LED count
        self.valid_combinations_by_count = {}
        self.configuration_counts_by_led = {}  # Number of configurations per LED count
        self._precompute_all_valid_combinations()
        
        if verbose:
            print(f"Initialized for grid {self.grid_size}x{self.grid_size} with {self.G.step} cm step size")
            print(f"LEDs: {led_min}-{led_max}")
            print(f"Heights: {height_min}-{height_max} (step: {height_step})")
            print(f"Number of possible heights: {len(self.possible_heights)}")
            print(f"Random seed: {seed if seed is not None else 'None (non-reproducible)'}")
            print(f"Sampling mode: {sampling_mode}")
            
            # Print configuration count
            self._print_configuration_per_led_count()
        
    def _print_configuration_per_led_count(self):
        """
        Prints the number of valid symmetric LED configurations for each possible LED count.

        :side effects:
            - Outputs the configuration count per LED count and the total number of configurations to stdout.

        :notes:
            - This method is intended for diagnostic and informational purposes, typically called during initialization if verbose mode is enabled.
            - The counts reflect precomputed, symmetry-constrained configurations for the current grid and LED parameters.
        """
        print(f"\nConfigurations per LED count:")
        total_configs = 0
        for led_count in sorted(self.configuration_counts_by_led.keys()):
            count = self.configuration_counts_by_led[led_count]
            total_configs += count
            print(f"  {led_count} LEDs: {count:,} configurations")
        print(f"Total configurations across all LED counts: {total_configs:,}\n")
        
    def _generate_height_range(self) -> List[float]:
        """
        Generates a list of all possible LED heights based on the minimum, maximum, and step size parameters.

        :returns:
            List[float]: A list of height values (in centimeters) that LEDs can be placed at, generated from
            ``height_min`` to ``height_max`` (inclusive) in increments of ``height_step``. Values are rounded to 6 decimal places.

        :notes:
            - Floating point precision is accounted for to ensure the maximum value is included if within tolerance.
            - The returned list is used for random sampling of LED heights in configuration generation.
        """
        heights = []
        current = self.height_min
        while current <= self.height_max + 1e-9:  # Account for floating point precision
            heights.append(round(current, 6))  # Round to avoid floating point issues
            current += self.height_step
        return heights
    
    def _classify_points(self):
        """
        Classifies all grid points into quadrants (Q1–Q4), axes (positive/negative x and y), and identifies the center.

        :side effects:
            - Populates the following instance attributes with lists of (x, y) tuples:
                - ``Q1_points``: Points in the first quadrant (x > 0, y > 0)
                - ``Q2_points``: Points in the second quadrant (x < 0, y > 0)
                - ``Q3_points``: Points in the third quadrant (x < 0, y < 0)
                - ``Q4_points``: Points in the fourth quadrant (x > 0, y < 0)
                - ``pos_x_axis``, ``neg_x_axis``: Points on the positive/negative x-axis (excluding center)
                - ``pos_y_axis``, ``neg_y_axis``: Points on the positive/negative y-axis (excluding center)
                - ``all_pos_axis``: All axis points with positive coordinates (excluding center)
            - The center point (0, 0) is handled separately and not included in any quadrant or axis list.

        :notes:
            - This classification is essential for enforcing symmetry and balance constraints in LED configuration generation.
            - Should be called during initialization before any configuration sampling.
        """
        for x in self.coords:
            for y in self.coords:
                if x == 0 and y == 0:
                    continue  # center handled separately
                elif x > 0 and y > 0:
                    self.Q1_points.append((x, y))
                elif x < 0 and y > 0:
                    self.Q2_points.append((x, y))
                elif x < 0 and y < 0:
                    self.Q3_points.append((x, y))
                elif x > 0 and y < 0:
                    self.Q4_points.append((x, y))
                elif y == 0 and x != 0:
                    if x > 0:
                        self.pos_x_axis.append((x, y))
                    else:
                        self.neg_x_axis.append((x, y))
                elif x == 0 and y != 0:
                    if y > 0:
                        self.pos_y_axis.append((x, y))
                    else:
                        self.neg_y_axis.append((x, y))
        self.all_pos_axis = self.pos_x_axis + self.pos_y_axis
        
    def _get_configuration_count(self, a: int, b: int, c: int) -> int:
        """
        Calculates the number of possible symmetric LED configurations for a specific combination of quadrant and axis pairs.

        :param a: Number of Type A pairs (Q1 <-> Q3), i.e., pairs of points symmetric about the center in quadrants 1 and 3. (int)
        :param b: Number of Type B pairs (Q2 <-> Q4), i.e., pairs of points symmetric about the center in quadrants 2 and 4. (int)
        :param c: Number of axis pairs, i.e., pairs of points symmetric about the center on the axes (positive/negative x or y, excluding center). (int)

        :returns:
            int: The total number of unique symmetric configurations that can be formed by selecting ``a`` Q1 pairs, ``b`` Q2 pairs, and ``c`` axis pairs.

        :notes:
            - The calculation is based on combinatorial selection of points from the classified quadrants and axes.
            - Used internally to determine the weight of each (a, b, c) combination for random sampling.
            - Assumes that the lists ``Q1_points``, ``Q2_points``, and ``all_pos_axis`` have been populated.
        """
        count = 1
        # Type A pairs: choose 'a' points from Q1
        if a > 0:
            q1_choices = len(list(combinations(self.Q1_points, a)))
            count *= q1_choices
        # Type B pairs: choose 'b' points from Q2  
        if b > 0:
            q2_choices = len(list(combinations(self.Q2_points, b)))
            count *= q2_choices
        # Axis pairs: choose 'c' points from axis points
        if c > 0:
            axis_choices = len(list(combinations(self.all_pos_axis, c)))
            count *= axis_choices
        return count
    
    def _precompute_all_valid_combinations(self):
        """
        Precomputes all valid symmetric LED configuration combinations for each possible LED count in the specified range.

        This method determines all valid (a, b, c) combinations, where:
        - a: Number of Type A pairs (Q1 <-> Q3)
        - b: Number of Type B pairs (Q2 <-> Q4)
        - c: Number of axis pairs (on axes, symmetric about center)
        for each LED count between ``self.led_min`` and ``self.led_max`` (inclusive), enforcing symmetry and quadrant balance.

        The method also calculates the number of configurations for each combination and stores the results for efficient sampling.

        :side effects:
            - Populates ``self.valid_combinations_by_count`` with valid combinations and their weights.
            - Populates ``self.configuration_counts_by_led`` with the number of configurations per LED count.
            - Sets ``self.total_configurations`` and ``self.led_count_weights`` for sampling.

        :notes:
            - This method is called during initialization and is essential for enabling efficient random sampling of valid configurations.
            - Assumes that the grid has already been classified into quadrants and axes.
        """
        total_configs = 0
        
        for led_count in range(self.led_min, self.led_max + 1):
            # use center for odd LED counts
            center_used = (led_count % 2 == 1)
            pair_count = led_count // 2
            if center_used:
                pair_count = (led_count - 1) // 2
                
            valid_combinations = []
            combo_config_counts = []
            target_weight = led_count / 4
            center_adjustment = 0.25 if center_used else 0.0
            
            for a in range(pair_count + 1):
                for b in range(pair_count + 1 - a):
                    c = pair_count - a - b
                    if c < 0:
                        continue
                    
                    q1_weight = a + 0.5 * c + center_adjustment
                    q2_weight = b + 0.5 * c + center_adjustment
                    
                    if (abs(q1_weight - target_weight) < 1e-6 and 
                        abs(q2_weight - target_weight) < 1e-6):
                        valid_combinations.append((a, b, c))
                        
                        # Calculate number of configurations for this combo
                        config_count = self._get_configuration_count(a, b, c)
                        combo_config_counts.append(config_count)
            
            if valid_combinations:
                total_for_led = sum(combo_config_counts)
                self.valid_combinations_by_count[led_count] = {
                    'combinations': valid_combinations,
                    'center_used': center_used,
                    'pair_count': pair_count,
                    'combo_weights': [count/total_for_led for count in combo_config_counts],
                    'total_configs': total_for_led
                }
                self.configuration_counts_by_led[led_count] = total_for_led
                total_configs += total_for_led
        
        # Store total for normalization
        self.total_configurations = total_configs
        
        # Calculate weights for LED count selection based on sampling mode
        self.led_count_weights = {}
        valid_led_counts = list(self.valid_combinations_by_count.keys())
        
        for led_count in valid_led_counts:
            if self.sampling_mode == "uniform":
                # Equal probability for each LED count
                self.led_count_weights[led_count] = 1.0 / len(valid_led_counts)
            else:
                # Weight by number of configurations
                self.led_count_weights[led_count] = (
                    self.configuration_counts_by_led[led_count] / total_configs
                )

    def get_quadrant_weight(self, point: Tuple[int, int]) -> Tuple[float, float, float, float]:
        """
        Determines the fractional contribution of a grid point to each of the four quadrants (Q1, Q2, Q3, Q4).

        :param point: Tuple[int, int]
            The (x, y) coordinates of the grid point.
        :returns:
            Tuple[float, float, float, float]:
                The weights for (Q1, Q2, Q3, Q4), each in [0, 1], summing to 1 for the center, or 1 for a single quadrant/axis.

        :notes:
            - Used for verifying quadrant balance in LED configurations.
            - Axis points contribute to two quadrants; the center contributes equally to all.
        """
        x, y = point
        if x == 0 and y == 0:
            return (0.25, 0.25, 0.25, 0.25)
        elif x > 0 and y > 0:  # Q1 point
            return (1.0, 0.0, 0.0, 0.0)
        elif x < 0 and y > 0:  # Q2 point  
            return (0.0, 1.0, 0.0, 0.0)
        elif x < 0 and y < 0:  # Q3 point
            return (0.0, 0.0, 1.0, 0.0)
        elif x > 0 and y < 0:  # Q4 point
            return (0.0, 0.0, 0.0, 1.0)
        elif y == 0 and x != 0:  # x-axis point
            if x > 0:  # positive x-axis
                return (0.5, 0.0, 0.0, 0.5)  # Contributes to Q1 and Q4
            else:  # negative x-axis  
                return (0.0, 0.5, 0.5, 0.0)  # FIXED: Contributes to Q2 and Q3
        elif x == 0 and y != 0:  # y-axis point
            if y > 0:  # positive y-axis
                return (0.5, 0.5, 0.0, 0.0)  # Contributes to Q1 and Q2
            else:  # negative y-axis
                return (0.0, 0.0, 0.5, 0.5)  # Contributes to Q3 and Q4
        return (0.0, 0.0, 0.0, 0.0)
    
    def calculate_quadrant_weights(self, configuration: List[Tuple[int, int]]) -> Tuple[float, float, float, float]:
        """
        Calculates the total fractional weight of LEDs in each quadrant for a given configuration.

        :param configuration: List[Tuple[int, int]]
            List of (x, y) grid points representing the LED configuration.
        :returns:
            Tuple[float, float, float, float]:
                The total weights for (Q1, Q2, Q3, Q4), summing to the total number of LEDs.

        :notes:
            - Used to verify that a configuration is balanced across all quadrants.
        """
        q1, q2, q3, q4 = 0.0, 0.0, 0.0, 0.0
        for point in configuration:
            w1, w2, w3, w4 = self.get_quadrant_weight(point)
            q1 += w1
            q2 += w2
            q3 += w3
            q4 += w4
        return q1, q2, q3, q4
    
    def verify_configuration(self, configuration: List[Tuple[int, int]]) -> bool:
        """
        Checks whether a given LED configuration satisfies all symmetry and balance constraints.

        :param configuration: List[Tuple[int, int]]
            List of (x, y) grid points representing the LED configuration.
        :returns:
            bool: True if the configuration is valid (correct LED count, point symmetry, and quadrant balance), False otherwise.

        :raises:
            None

        :notes:
            - The configuration must have a valid LED count, be symmetric about the center, and be balanced across all quadrants.
        """
        led_count = len(configuration)
        is_valid = False

        # Check if LED count is in valid range
        if not (self.led_min <= led_count <= self.led_max):
            return is_valid
        
        # Check point symmetry
        config_set = set(configuration)
        for point in configuration:
            symmetric_point = (-point[0], -point[1])
            if symmetric_point not in config_set:
                return is_valid
        
        # Check quadrant balance
        q1, q2, q3, q4 = self.calculate_quadrant_weights(configuration)
        target_weight = led_count / 4

        is_valid = (abs(q1 - target_weight) < 1e-6 and 
                abs(q2 - target_weight) < 1e-6 and
                abs(q3 - target_weight) < 1e-6 and
                abs(q4 - target_weight) < 1e-6)
        return is_valid
    
    def set_seed(self, seed: int):
        """
        Sets the random seed for reproducible sampling.

        :param seed: int
            The random seed value to use for both Python's random module and NumPy.

        :side effects:
            - Updates the internal random state and NumPy's global seed.
            - Prints the seed value to stdout.

        :notes:
            - Use this method to ensure reproducibility of generated configurations.
        """
        self.seed = seed
        self._random_state = random.Random(seed)
        if seed is not None:
            np.random.seed(seed)
        print(f"Random seed set to: {seed}")
    
    def _random_choice(self, seq, weights=None):
        """
        Selects a single random element from a sequence, optionally using weights, with the sampler's random state.

        :param seq: Sequence
            The sequence to choose from.
        :param weights: Optional[List[float]]
            Optional weights for each element in the sequence.
        :returns:
            Any: The randomly selected element from the sequence.

        :notes:
            - Ensures reproducibility by using the sampler's internal random state.
        """
        if weights is None:
            return self._random_state.choice(seq)
        else:
            return self._random_state.choices(seq, weights=weights, k=1)[0]
    
    def _random_choices(self, population, weights=None, k=1):
        """
        Selects multiple random elements from a population, optionally using weights, with the sampler's random state.

        :param population: Sequence
            The population to choose from.
        :param weights: Optional[List[float]]
            Optional weights for each element in the population.
        :param k: int
            Number of elements to select (default: 1).
        :returns:
            List[Any]: List of randomly selected elements.

        :notes:
            - Ensures reproducibility by using the sampler's internal random state.
        """
        if weights is None:
            return self._random_state.choices(population, k=k)
        else:
            return self._random_state.choices(population, weights=weights, k=k)
    
    def _random_sample(self, population, k):
        """
        Draws a random sample of unique elements from a population using the sampler's random state.

        :param population: Sequence
            The population to sample from.
        :param k: int
            Number of unique samples to draw.
        :returns:
            List[Any]: List of sampled elements.

        :raises ValueError:
            If k is larger than the population size.

        :notes:
            - Ensures reproducibility by using the sampler's internal random state.
        """
        return self._random_state.sample(population, k)
    
    def generate_random_configuration(self, 
                                target_led_count: Optional[int] = None,
                                target_height: Optional[float] = None) -> Optional[Dict[str, Any]]:
        """
        Generates a single random valid LED configuration and associated height, enforcing symmetry and quadrant balance.

        :param target_led_count: Optional[int]
            If specified, generates a configuration with this exact number of LEDs. If None, the LED count is sampled according to the sampler's distribution.
        :param target_height: Optional[float]
            If specified, selects the LED height closest to this value. If None, height is sampled from possible values.
        :returns:
            Optional[Dict[str, Any]]: A dictionary containing configuration details (positions, LED count, height, symmetry tuple, and seed), or None if generation failed after maximum attempts.

        :side effects:
            - May print debug information to stdout if configuration generation fails or is invalid.

        :notes:
            - Attempts up to 1000 times to generate a valid configuration.
            - Uses the sampler's random state for reproducibility.
            - The returned dictionary includes the configuration, LED count, height, symmetry tuple (a, b, c), and seed used.
        """
        max_attempts = 1000
        
        # Determine LED count
        if target_led_count is not None:
            if target_led_count not in self.valid_combinations_by_count:
                print(f"DEBUG: Target LED {target_led_count} not in valid combinations")
                return None
            valid_led_counts = [target_led_count]
            led_weights = None
        else:
            valid_led_counts = list(self.valid_combinations_by_count.keys())
            led_weights = [self.led_count_weights[lc] for lc in valid_led_counts]
        
        # Determine height
        if target_height is not None:
            height = min(self.possible_heights, key=lambda h: abs(h - target_height))
        else:
            height = self._random_choice(self.possible_heights)
        
        led_count = None
        for attempt in range(max_attempts):
            try:
                # Randomly choose LED count
                if target_led_count is None:
                    led_count = self._random_choice(valid_led_counts, weights=led_weights)
                else:
                    led_count = target_led_count
                
                combo_info = self.valid_combinations_by_count[led_count]
                center_used = combo_info['center_used']
                pair_count = combo_info['pair_count']
                
                # Randomly choose a valid (a,b,c) combination with appropriate weights
                a, b, c = self._random_choice(
                    combo_info['combinations'], 
                    weights=combo_info['combo_weights']
                )
                
                config = set()
                
                # Add center if needed
                if center_used:
                    config.add(self.center)
                
                # Add Type A pairs (Q1 -> Q3)
                if a > 0:
                    q1_selection = self._random_sample(self.Q1_points, a)
                    for pt in q1_selection:
                        config.add(pt)
                        config.add((-pt[0], -pt[1]))
                
                # Add Type B pairs (Q2 -> Q4)
                if b > 0:
                    q2_selection = self._random_sample(self.Q2_points, b)
                    for pt in q2_selection:
                        config.add(pt)
                        config.add((-pt[0], -pt[1]))
                
                # Add axis pairs
                if c > 0:
                    axis_selection = self._random_sample(self.all_pos_axis, c)
                    for axis_point in axis_selection:
                        config.add(axis_point)
                        config.add((-axis_point[0], -axis_point[1]))
                
                # Verify we have the correct number of points
                if len(config) != led_count:
                    if attempt == 0:  # Only print once to avoid spam
                        print(f"DEBUG: LED {led_count} - Wrong point count: got {len(config)}, expected {led_count}")
                    continue
                
                sorted_config = sorted(config)
                if self.verify_configuration(sorted_config):
                    return {
                        'configuration': sorted_config,
                        'led_count': led_count,
                        'height': height,
                        'a_b_c': (a, b, c),
                        # 'positions': sorted_config,
                        'seed_used': self.seed
                    }
                else:
                    if attempt == 0:  # Only print once to avoid spam
                        print(f"DEBUG: LED {led_count} - Failed verification")
                        q1, q2, q3, q4 = self.calculate_quadrant_weights(sorted_config)
                        target = led_count / 4
                        print(f"  Weights: Q1={q1}, Q2={q2}, Q3={q3}, Q4={q4}, target={target}")
            except (ValueError, KeyError) as e:
                if attempt == 0:  # Only print once to avoid spam
                    print(f"DEBUG: LED {led_count if led_count is not None else 'Unknown'} - Exception: {e}")
                continue
        print(f"DEBUG: LED {led_count if led_count is not None else 'Unknown'} - Failed after {max_attempts} attempts")
        return None

    def generate_samples(self, num_samples: int, 
                        led_distribution: Optional[Dict[int, float]] = None,
                        height_distribution: Optional[Dict[float, float]] = None,
                        use_parallel: bool = True, 
                        max_workers: int = None,
                        worker_seeds: Optional[List[int]] = None,
                        target_samples_per_led: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Generates multiple random LED configurations with variable LED counts and heights, supporting parallel processing and custom distributions.

        This method creates a specified number of random, symmetric LED configurations on a square grid, optionally using user-defined distributions for LED counts and heights. Sampling can be performed in parallel for efficiency. If `target_samples_per_led` is specified, the method generates approximately that many samples for each valid LED count.

        :param num_samples: Number of samples to generate. (int)
        :param led_distribution: Optional; probability distribution for LED counts, as a dictionary mapping LED count (int) to probability (float). If None, uses the sampler's default distribution. (Optional[Dict[int, float]])
        :param height_distribution: Optional; probability distribution for LED heights, as a dictionary mapping height (float) to probability (float). If None, uses uniform sampling over allowed heights. (Optional[Dict[float, float]])
        :param use_parallel: Whether to use parallel processing for sample generation. (bool, default: True)
        :param max_workers: Maximum number of parallel worker processes. If None, uses the number of CPU cores. (Optional[int])
        :param worker_seeds: Optional list of random seeds for each worker to ensure reproducibility in parallel mode. (Optional[List[int]])
        :param target_samples_per_led: If specified, generates approximately this many samples for each valid LED count, overriding `num_samples`. (Optional[int])

        :returns: List of configuration dictionaries, each containing LED positions, count, height, and other metadata. (List[Dict[str, Any]])

        :raises ValueError: If parameter constraints are violated or grid is not properly initialized.

        :side effects:
            - May print progress and distribution information to stdout if verbose mode is enabled.
            - May use parallel processing and spawn worker processes.

        :notes:
            - If `target_samples_per_led` is set, the total number of samples is determined by the number of valid LED counts times this value.
            - Parallel processing is recommended for large sample counts.
            - All random sampling is reproducible if a seed is provided.
        """
        # If target_samples_per_led is specified, override the sampling strategy
        if target_samples_per_led is not None:
            return self._generate_samples_fixed_per_led(target_samples_per_led, height_distribution, use_parallel)
        
        # Print initialization details
        if self.verbose:
            print(f"Generating {num_samples} Monte Carlo samples with:")
            print(f"  LEDs: {self.led_min}-{self.led_max}")
            print(f"  Heights: {self.height_min}-{self.height_max} (step: {self.height_step})")
            print(f"  Seed: {self.seed if self.seed is not None else 'None'}")
            print(f"  Sampling mode: {self.sampling_mode}")
            
            if led_distribution:
                print(f"  LED distribution: {led_distribution}")
            if height_distribution:
                print(f"  Height distribution: {height_distribution}")
        
        if use_parallel:
            return self._generate_samples_parallel(num_samples, led_distribution, 
                                                 height_distribution, max_workers, worker_seeds)
        else:
            return self._generate_samples_sequential(num_samples, led_distribution, height_distribution)

    def _generate_samples_fixed_per_led(self, samples_per_led: int,
                                      height_distribution: Optional[Dict[float, float]] = None,
                                      use_parallel: bool = True) -> List[Dict[str, Any]]:
        """
        Generates approximately the same number of random LED configurations for each valid LED count.

        This method ensures that each allowed LED count is represented by roughly `samples_per_led` samples, optionally using a custom height distribution and parallel processing.

        :param samples_per_led: Number of samples to generate for each valid LED count. (int)
        :param height_distribution: Optional; probability distribution for LED heights, as a dictionary mapping height (float) to probability (float). If None, uses uniform sampling. (Optional[Dict[float, float]])
        :param use_parallel: Whether to use parallel processing for sample generation. (bool, default: True)

        :returns: List of configuration dictionaries, each containing LED positions, count, height, and other metadata. (List[Dict[str, Any]])

        :side effects:
            - Prints progress and distribution information to stdout if verbose mode is enabled.
            - May use parallel processing and spawn worker processes.

        :notes:
            - The total number of samples generated is `samples_per_led` times the number of valid LED counts.
            - Each LED count is sampled independently.
        """
        valid_led_counts = sorted(self.valid_combinations_by_count.keys())
        total_samples = samples_per_led * len(valid_led_counts)
        
        if self.verbose:
            print(f"Generating samples with fixed {samples_per_led} samples per LED count") 
            print(f"Total target samples: {total_samples}")
            print(f"LED counts: {valid_led_counts}")
        
        all_samples = []
        
        for led_count in valid_led_counts:
            if self.verbose:
                print(f"\nGenerating {samples_per_led} samples for {led_count} LEDs...")
            
            if use_parallel:
                led_samples = self._generate_samples_parallel(
                    samples_per_led, 
                    led_distribution={led_count: 1.0},  # Force this LED count
                    height_distribution=height_distribution,
                    max_workers=min(os.cpu_count(), samples_per_led)
                )
            else:
                led_samples = self._generate_samples_sequential(
                    samples_per_led,
                    led_distribution={led_count: 1.0},  # Force this LED count
                    height_distribution=height_distribution
                )
            
            all_samples.extend(led_samples)
            if self.verbose:
                print(f"  Generated {len(led_samples)} valid samples for {led_count} LEDs")
        
        # Analyze final distribution
        if self.verbose:
            led_count_dist = defaultdict(int)
            for sample in all_samples:
                led_count_dist[sample['led_count']] += 1

            print(f"\nFinal LED count distribution:")
            for led_count in sorted(led_count_dist.keys()):
                count = led_count_dist[led_count]
                percentage = (count / len(all_samples)) * 100
                print(f"  {led_count} LEDs: {count} samples ({percentage:.1f}%)")
        
        return all_samples

    def _generate_samples_sequential(self, num_samples: int, 
                                   led_distribution: Optional[Dict[int, float]] = None,
                                   height_distribution: Optional[Dict[float, float]] = None) -> List[Dict[str, Any]]:
        """
        Generates random LED configurations sequentially (single-threaded), optionally using custom LED count and height distributions.

        :param num_samples: Number of samples to generate. (int)
        :param led_distribution: Optional; probability distribution for LED counts, as a dictionary mapping LED count (int) to probability (float). If None, uses the sampler's default distribution. (Optional[Dict[int, float]])
        :param height_distribution: Optional; probability distribution for LED heights, as a dictionary mapping height (float) to probability (float). If None, uses uniform sampling. (Optional[Dict[float, float]])

        :returns: List of configuration dictionaries, each containing LED positions, count, height, and other metadata. (List[Dict[str, Any]])

        :side effects:
            - Prints progress and distribution information to stdout if verbose mode is enabled.

        :notes:
            - This method is suitable for small to moderate sample counts or debugging.
            - All random sampling is reproducible if a seed is provided.
        """
        samples = []
        start_time = time.time()
        led_count_dist = defaultdict(int)
        height_dist = defaultdict(int)
        
        for i in range(num_samples):
            # Determine target LED count based on distribution
            target_led_count = None
            if led_distribution:
                target_led_count = self._random_choices(
                    list(led_distribution.keys()),
                    weights=list(led_distribution.values())
                )[0]
            
            # Determine target height based on distribution
            target_height = None
            if height_distribution:
                target_height = self._random_choices(
                    list(height_distribution.keys()),
                    weights=list(height_distribution.values())
                )[0]
            
            config_data = self.generate_random_configuration(target_led_count, target_height)
            if config_data is not None:
                samples.append(config_data)
                led_count_dist[config_data['led_count']] += 1
                height_dist[config_data['height']] += 1
        
        # Print distributions
        self._print_distributions(led_count_dist, height_dist, num_samples)
        
        return samples

    def _generate_samples_parallel(self, num_samples: int,
                                 led_distribution: Optional[Dict[int, float]] = None,
                                 height_distribution: Optional[Dict[float, float]] = None,
                                 max_workers: int = None,
                                 worker_seeds: Optional[List[int]] = None) -> List[Dict[str, Any]]:
        """
        Generates random LED configurations in parallel using multiple worker processes, optionally with custom LED count and height distributions.

        :param num_samples: Number of samples to generate. (int)
        :param led_distribution: Optional; probability distribution for LED counts, as a dictionary mapping LED count (int) to probability (float). If None, uses the sampler's default distribution. (Optional[Dict[int, float]])
        :param height_distribution: Optional; probability distribution for LED heights, as a dictionary mapping height (float) to probability (float). If None, uses uniform sampling. (Optional[Dict[float, float]])
        :param max_workers: Maximum number of parallel worker processes. If None, uses the number of CPU cores. (Optional[int])
        :param worker_seeds: Optional list of random seeds for each worker to ensure reproducibility. (Optional[List[int]])

        :returns: List of configuration dictionaries, each containing LED positions, count, height, and other metadata. (List[Dict[str, Any]])

        :side effects:
            - Prints progress and distribution information to stdout.
            - Spawns worker processes for parallel execution.

        :notes:
            - Parallel processing is recommended for large sample counts.
            - Each worker can be seeded for reproducibility.
        """
        if max_workers is None:
            max_workers = os.cpu_count()
        
        # Generate worker seeds if not provided
        if worker_seeds is None and self.seed is not None:
            # Generate deterministic seeds for workers based on main seed
            worker_seeds = [self.seed + i for i in range(max_workers)]
        
        print(f"Using {max_workers} parallel workers")
        if worker_seeds:
            print(f"Worker seeds: {worker_seeds}")
        
        # Split work among workers
        samples_per_worker = [num_samples // max_workers] * max_workers
        for i in range(num_samples % max_workers):
            samples_per_worker[i] += 1
        
        start_time = time.time()
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit tasks to workers
            futures = []
            for i, worker_samples in enumerate(samples_per_worker):
                if worker_samples > 0:
                    worker_seed = worker_seeds[i] if worker_seeds else None
                    future = executor.submit(self._worker_generate_samples, worker_samples, 
                                           led_distribution, height_distribution, worker_seed)
                    futures.append(future)
            
            # Collect results
            all_samples = []
            led_count_dist = defaultdict(int)
            height_dist = defaultdict(int)
            completed = 0
            
            for future in concurrent.futures.as_completed(futures):
                worker_samples = future.result()
                all_samples.extend(worker_samples)
                
                # Update distribution counts
                for sample in worker_samples:
                    led_count_dist[sample['led_count']] += 1
                    height_dist[sample['height']] += 1
                
                completed += 1
                elapsed = time.time() - start_time
                total_generated = len(all_samples)
                rate = total_generated / elapsed if elapsed > 0 else 0
                print(f"Worker {completed}/{len(futures)} completed: {len(worker_samples)} samples "
                      f"(Total: {total_generated}, {rate:.1f} samples/sec)")
        
        # Print distributions
        self._print_distributions(led_count_dist, height_dist, num_samples)
        
        return all_samples

    def _worker_generate_samples(self, num_samples: int,
                               led_distribution: Optional[Dict[int, float]] = None,
                               height_distribution: Optional[Dict[float, float]] = None,
                               worker_seed: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Generates random LED configurations in a worker process for parallel sample generation.

        This function is intended to be called by a parallel worker and uses its own random state for reproducibility.

        :param num_samples: Number of samples to generate in this worker. (int)
        :param led_distribution: Optional; probability distribution for LED counts, as a dictionary mapping LED count (int) to probability (float). If None, uses the sampler's default distribution. (Optional[Dict[int, float]])
        :param height_distribution: Optional; probability distribution for LED heights, as a dictionary mapping height (float) to probability (float). If None, uses uniform sampling. (Optional[Dict[float, float]])
        :param worker_seed: Random seed for this worker's random state. (Optional[int])

        :returns: List of configuration dictionaries, each containing LED positions, count, height, and other metadata. (List[Dict[str, Any]])

        :side effects:
            - Uses a worker-specific random state for reproducibility.

        :notes:
            - This method is not intended to be called directly by users; it is used internally for parallel sample generation.
        """
        # Create a new random state for this worker
        worker_random = random.Random(worker_seed)
        
        def worker_random_choice(seq, weights=None):
            if weights is None:
                return worker_random.choice(seq)
            else:
                return worker_random.choices(seq, weights=weights, k=1)[0]
        
        def worker_random_choices(population, weights=None, k=1):
            if weights is None:
                return worker_random.choices(population, k=k)
            else:
                return worker_random.choices(population, weights=weights, k=k)
        
        def worker_random_sample(population, k):
            return worker_random.sample(population, k)
        
        worker_samples = []
        for _ in range(num_samples):
            # Determine target LED count based on distribution
            target_led_count = None
            if led_distribution:
                target_led_count = worker_random_choices(
                    list(led_distribution.keys()),
                    weights=list(led_distribution.values())
                )[0]
            
            # Determine target height based on distribution
            target_height = None
            if height_distribution:
                target_height = worker_random_choices(
                    list(height_distribution.keys()),
                    weights=list(height_distribution.values())
                )[0]
            
            # For the actual configuration generation, we'll use the original method
            # but with worker-specific random functions
            config_data = self._generate_single_with_custom_random(
                target_led_count, target_height, worker_random_choice, worker_random_sample
            )
            if config_data is not None:
                config_data['worker_seed'] = worker_seed
                config_data['configuration'] = self._get_actual_coordinates(config_data['configuration'])
                worker_samples.append(config_data)
        
        return worker_samples

    def _generate_single_with_custom_random(self, target_led_count, target_height,
                                          random_choice_fn, random_sample_fn):
        """
        Generates a single random LED configuration using custom random functions, for use in parallel workers.

        This method is used internally to allow parallel workers to generate configurations with their own random state.

        :param target_led_count: Target number of LEDs for the configuration. If None, sampled according to distribution. (int or None)
        :param target_height: Target LED height. If None, sampled from possible heights. (float or None)
        :param random_choice_fn: Function for random selection from a sequence, supporting optional weights. (Callable)
        :param random_sample_fn: Function for random sampling of unique elements from a sequence. (Callable)

        :returns: Configuration dictionary with LED positions, count, height, and metadata, or None if generation fails. (Dict[str, Any] or None)

        :notes:
            - Used only in parallel sample generation; not intended for direct user calls.
            - Attempts up to 1000 times to generate a valid configuration.
        """
        max_attempts = 1000
        
        # Determine LED count
        if target_led_count is not None:
            if target_led_count not in self.valid_combinations_by_count:
                return None
            valid_led_counts = [target_led_count]
            led_weights = None
        else:
            valid_led_counts = list(self.valid_combinations_by_count.keys())
            led_weights = [self.led_count_weights[lc] for lc in valid_led_counts]
        
        # Determine height
        if target_height is not None:
            height = min(self.possible_heights, key=lambda h: abs(h - target_height))
        else:
            height = random_choice_fn(self.possible_heights)
        
        for attempt in range(max_attempts):
            try:
                # Randomly choose LED count
                if target_led_count is None:
                    led_count = random_choice_fn(valid_led_counts, weights=led_weights)
                else:
                    led_count = target_led_count
                    
                combo_info = self.valid_combinations_by_count[led_count]
                center_used = combo_info['center_used']
                pair_count = combo_info['pair_count']
                
                # Randomly choose a valid (a,b,c) combination with appropriate weights
                a, b, c = random_choice_fn(
                    combo_info['combinations'], 
                    weights=combo_info['combo_weights']
                )
                
                config = set()
                
                if center_used:
                    config.add(self.center)
                
                if a > 0:
                    q1_selection = random_sample_fn(self.Q1_points, a)
                    for pt in q1_selection:
                        config.add(pt)
                        config.add((-pt[0], -pt[1]))
                
                if b > 0:
                    q2_selection = random_sample_fn(self.Q2_points, b)
                    for pt in q2_selection:
                        config.add(pt)
                        config.add((-pt[0], -pt[1]))
                
                if c > 0:
                    axis_selection = random_sample_fn(self.all_pos_axis, c)
                    for axis_point in axis_selection:
                        config.add(axis_point)
                        config.add((-axis_point[0], -axis_point[1]))
                
                if len(config) == led_count:
                    sorted_config = sorted(config)
                    if self.verify_configuration(sorted_config):
                        return {
                            'configuration': sorted_config,
                            'led_count': led_count,
                            'height': height,
                            'a_b_c': (a, b, c),
                            # 'positions': sorted_config
                        }
                        
            except (ValueError, KeyError):
                continue
        
        return None
    
    def _get_actual_coordinates(self, samples: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Converts grid-based LED coordinates in sample dictionaries to real-world coordinates using the grid step size.

        :param samples: List of sample dictionaries, each containing a 'configuration' key with grid coordinates. (List[Dict[str, Any]])

        :returns: List of sample dictionaries with 'configuration' values converted to real-world coordinates. (List[Dict[str, Any]])

        :notes:
            - This method modifies the input dictionaries in-place.
            - The conversion multiplies each coordinate by the grid step size.
        """
        for sample in samples:
            sample['configuration'] = [(x * self.grid_step, y * self.grid_step) for (x, y) in sample['configuration']]
        return samples

    def _print_distributions(self, led_count_dist: Dict[int, int], 
                           height_dist: Dict[float, int], total_samples: int):
        """
        Prints the distribution of LED counts and heights among generated samples.

        :param led_count_dist: Dictionary mapping LED count (int) to number of samples. (Dict[int, int])
        :param height_dist: Dictionary mapping height (float) to number of samples. (Dict[float, int])
        :param total_samples: Total number of samples generated. (int)

        :side effects:
            - Outputs distribution statistics to stdout.

        :notes:
            - Intended for diagnostic and informational purposes.
        """
        if not led_count_dist:
            return
            
        print(f"\nLED Count Distribution:")
        for led_count in sorted(led_count_dist.keys()):
            count = led_count_dist[led_count]
            percentage = (count / total_samples) * 100
            print(f"  {led_count} LEDs: {count} samples ({percentage:.1f}%)")
        
        print(f"\nHeight Distribution:")
        for height in sorted(height_dist.keys()):
            count = height_dist[height]
            percentage = (count / total_samples) * 100
            print(f"  Height {height}: {count} samples ({percentage:.1f}%)")

    def generate_unique_samples(self, num_samples: int = 100, 
                              max_attempts_factor: int = 10,
                              led_distribution: Optional[Dict[int, float]] = None,
                              height_distribution: Optional[Dict[float, float]] = None,
                              use_parallel: bool = True,
                              target_samples_per_led: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Generates a specified number of unique, random, symmetric LED configurations with variable LED counts and heights.

        This method attempts to generate ``num_samples`` unique configurations, each satisfying symmetry and quadrant balance constraints, with LED counts and heights optionally sampled from user-provided distributions. Uniqueness is enforced on the combination of LED positions and heights. Sampling can be performed in parallel for efficiency. If ``target_samples_per_led`` is specified, the method generates approximately that many unique samples for each valid LED count.

        :param num_samples: Number of unique samples to generate. (int, default: 100)
        :param max_attempts_factor: Multiplier for the maximum number of attempts allowed (max_attempts = num_samples * max_attempts_factor). (int, default: 10)
        :param led_distribution: Optional; probability distribution for LED counts, as a dictionary mapping LED count (int) to probability (float). If None, uses the sampler's default distribution. (Optional[Dict[int, float]], default: None)
        :param height_distribution: Optional; probability distribution for LED heights, as a dictionary mapping height (float) to probability (float). If None, uses uniform sampling over allowed heights. (Optional[Dict[float, float]], default: None)
        :param use_parallel: Whether to use parallel processing for sample generation. (bool, default: True)
        :param target_samples_per_led: If specified, generates approximately this many unique samples for each valid LED count, overriding ``num_samples``. (Optional[int], default: None)

        :returns:
            List[Dict[str, Any]]: A list of unique configuration dictionaries, each containing LED positions, count, height, and metadata.

        :side effects:
            - Prints progress and distribution information to stdout.
            - May use parallel processing and spawn worker processes.

        :notes:
            - Uniqueness is enforced on the combination of LED positions and heights.
            - If ``target_samples_per_led`` is set, the total number of samples is determined by the number of valid LED counts times this value.
            - The method will stop after ``num_samples * max_attempts_factor`` attempts if not enough unique samples are found.
            - All random sampling is reproducible if a seed is provided.
        """
        # If target_samples_per_led is specified, use the fixed-per-LED approach
        if target_samples_per_led is not None:
            return self._generate_unique_samples_fixed_per_led(target_samples_per_led, height_distribution, use_parallel, max_attempts_factor)
            
        print(f"Generating {num_samples} UNIQUE Monte Carlo samples with:")
        print(f"  LEDs: {self.led_min}-{self.led_max}")
        print(f"  Heights: {self.height_min}-{self.height_max}")
        print(f"  Seed: {self.seed if self.seed is not None else 'None'}")
        print(f"  Sampling mode: {self.sampling_mode}")
        
        unique_samples = set()
        all_samples = []
        max_attempts = num_samples * max_attempts_factor
        attempts = 0
        
        start_time = time.time()
        led_count_dist = defaultdict(int)
        height_dist = defaultdict(int)
        
        while len(unique_samples) < num_samples and attempts < max_attempts:
            # Determine target LED count based on distribution
            target_led_count = None
            if led_distribution:
                target_led_count = self._random_choices(
                    list(led_distribution.keys()),
                    weights=list(led_distribution.values())
                )[0]
            
            # Determine target height based on distribution
            target_height = None
            if height_distribution:
                target_height = self._random_choices(
                    list(height_distribution.keys()),
                    weights=list(height_distribution.values())
                )[0]
            
            config_data = self.generate_random_configuration(target_led_count, target_height)
            attempts += 1
            
            if config_data is not None:
                # Create unique key based on positions AND height
                config_key = (tuple(tuple(point) for point in config_data['configuration']), 
                            config_data['height'])
                
                if config_key not in unique_samples:
                    unique_samples.add(config_key)
                    all_samples.append(config_data)
                    led_count_dist[config_data['led_count']] += 1
                    height_dist[config_data['height']] += 1
        
        # Print distributions
        self._print_distributions(led_count_dist, height_dist, len(all_samples))
        
        print(f"Generated {len(all_samples)} unique samples from {attempts} attempts "
              f"({len(all_samples)/attempts*100:.1f}% success rate)")
        
        # multiply the sample positions by grid step size
        for sample in all_samples:
            sample['configuration'] = [(x * self.grid_step, y * self.grid_step) for (x, y) in sample['configuration']]
        
        return all_samples

    def _generate_unique_samples_fixed_per_led(self, samples_per_led: int,
                                             height_distribution: Optional[Dict[float, float]] = None,
                                             use_parallel: bool = True,
                                             max_attempts_factor: int = 10) -> List[Dict[str, Any]]:
        """
        Generates approximately the same number of unique, random, symmetric LED configurations for each valid LED count.

        This method ensures that each allowed LED count is represented by roughly ``samples_per_led`` unique samples, optionally using a custom height distribution and parallel processing. Uniqueness is enforced on the combination of LED positions and heights.

        :param samples_per_led: Number of unique samples to generate per valid LED count. (int)
        :param height_distribution: Optional; probability distribution for LED heights, as a dictionary mapping height (float) to probability (float). If None, uses uniform sampling. (Optional[Dict[float, float]], default: None)
        :param use_parallel: Whether to use parallel processing for sample generation. (bool, default: True)
        :param max_attempts_factor: Multiplier for the maximum number of attempts allowed per LED count (max_attempts = samples_per_led * max_attempts_factor). (int, default: 10)

        :returns:
            List[Dict[str, Any]]: A list of unique configuration dictionaries, each containing LED positions, count, height, and metadata.

        :side effects:
            - Prints progress and distribution information to stdout.
            - May use parallel processing and spawn worker processes.

        :notes:
            - The total number of samples generated is ``samples_per_led`` times the number of valid LED counts.
            - Each LED count is sampled independently.
            - Uniqueness is enforced on the combination of LED positions and heights.
        """
        valid_led_counts = sorted(self.valid_combinations_by_count.keys())
        total_target_samples = samples_per_led * len(valid_led_counts)
        
        print(f"Generating {samples_per_led} unique samples for each of {len(valid_led_counts)} LED counts")
        print(f"Total target unique samples: {total_target_samples}")
        print(f"LED counts: {valid_led_counts}")
        
        all_samples = []
        
        for led_count in valid_led_counts:
            print(f"\nGenerating {samples_per_led} unique samples for {led_count} LEDs...")
            
            led_samples = self.generate_unique_samples(
                samples_per_led,
                max_attempts_factor=max_attempts_factor,
                led_distribution={led_count: 1.0},  # Force this LED count
                height_distribution=height_distribution,
                use_parallel=use_parallel
            )
            
            all_samples.extend(led_samples)
            print(f"Generated {len(led_samples)}/{samples_per_led} unique samples for {led_count} LEDs")
        
        # Analyze final distribution
        led_count_dist = defaultdict(int)
        for sample in all_samples:
            led_count_dist[sample['led_count']] += 1
        
        print(f"\nFinal LED Count Distribution:")
        for led_count in sorted(led_count_dist.keys()):
            count = led_count_dist[led_count]
            percentage = (count / len(all_samples)) * 100
            print(f"  {led_count} LEDs: {count} samples ({percentage:.1f}%)")
        
        return all_samples


def save_configurations(configurations: List[Dict[str, Any]], 
                       sampler: MonteCarloSampler,
                       filename: str):
    """
    Saves a list of LED configuration dictionaries to a compressed NumPy NPZ file, including metadata such as heights, LED counts, and generation parameters.

    The function extracts relevant data from the provided configurations and stores them in a compressed NPZ file for efficient storage and later retrieval. Metadata about the grid, LED range, height range, and random seed are also included for reproducibility and documentation.

    :param configurations: List of configuration dictionaries, each containing at least 'configuration', 'led_count', and 'height'. (List[Dict[str, Any]])
    :param sampler: The MonteCarloSampler instance used to generate the configurations. (MonteCarloSampler)
    :param filename: Output filename for the compressed NPZ file. (str)

    :returns:
        None

    :side effects:
        - Writes a compressed NPZ file to disk at the specified filename.
        - Prints a summary message to stdout upon successful save.

    :notes:
        - The NPZ file contains arrays for configurations, LED counts, heights, and metadata fields for reproducibility.
        - The function assumes all configurations are compatible with the provided sampler's grid and parameter ranges.
    """
    # Extract data for storage
    n_configs = len(configurations)
    
    # Extract simple data types using fromiter (fastest for numeric data)
    led_counts = np.fromiter((sample['led_count'] for sample in configurations), 
                           dtype=np.int32, count=n_configs)
    height = np.fromiter((sample['height'] for sample in configurations), 
                        dtype=np.float64, count=n_configs)
    configs = np.fromiter((sample['configuration'] for sample in configurations), 
                        dtype=object, count=n_configs)
    
    np.savez_compressed(
        filename,
        configurations=configs,
        led_counts=led_counts,
        height=height,
        total_count=len(configurations),
        grid_size=sampler.grid_size,
        led_min=sampler.led_min,
        led_max=sampler.led_max,
        height_min=sampler.height_min,
        height_max=sampler.height_max,
        seed=sampler.seed if sampler.seed is not None else -1,
        grid_step=sampler.grid_step,
        sampler="MonteCarloSampler",
        description=f"Variable LEDs {sampler.led_min}-{sampler.led_max} with heights {sampler.height_min}-{sampler.height_max} on {sampler.grid_size}x{sampler.grid_size} grid (seed: {sampler.seed}). Generated by VariableLEDHeightMonteCarloGenerator."
    )
    
    print(f"Saved {len(configurations)} configurations to {filename}")

def analyze_configurations(configurations: List[Dict[str, Any]], 
                          generator: MonteCarloSampler):
    """
    Analyzes a list of generated LED configurations and prints detailed statistics, including distributions, uniqueness, and validity checks.

    This function provides a summary of the generated configurations, including the distribution of LED counts and heights, the number of unique configurations (by position and by position+height), and a validity check on a sample of configurations using the generator's verification method.

    :param configurations: List of configuration dictionaries, each containing at least 'configuration', 'led_count', and 'height'. (List[Dict[str, Any]])
    :param generator: The MonteCarloSampler instance used to generate or verify the configurations. (MonteCarloSampler)

    :returns:
        None

    :side effects:
        - Prints analysis results and statistics to stdout.

    :notes:
        - Uniqueness is checked both for (positions + height) and for positions only.
        - Validity is checked for up to 10 randomly selected configurations using the generator's verification method.
        - The function assumes all configurations are compatible with the provided generator's grid and parameter ranges.
    """
    if not configurations:
        print("No configurations to analyze")
        return
    
    print(f"\n=== ANALYSIS ===")
    print(f"Total configurations: {len(configurations)}")
    
    # LED count distribution
    led_count_dist = defaultdict(int)
    height_dist = defaultdict(int)
    for config_data in configurations:
        led_count_dist[config_data['led_count']] += 1
        height_dist[config_data['height']] += 1
    
    print(f"LED Count Distribution:")
    for led_count in sorted(led_count_dist.keys()):
        count = led_count_dist[led_count]
        percentage = (count / len(configurations)) * 100
        print(f"  {led_count} LEDs: {count} configurations ({percentage:.1f}%)")
    
    print(f"\nHeight Distribution:")
    for height in sorted(height_dist.keys()):
        count = height_dist[height]
        percentage = (count / len(configurations)) * 100
        print(f"  Height {height}: {count} configurations ({percentage:.1f}%)")
    
    # Check uniqueness (positions + height)
    config_keys = [(tuple(tuple(point) for point in config_data['configuration']), 
                   config_data['height']) for config_data in configurations]
    unique_count = len(set(config_keys))
    print(f"Unique configurations (position + height): {unique_count}")
    
    # Check uniqueness (positions only)
    position_keys = [tuple(tuple(point) for point in config_data['configuration']) 
                    for config_data in configurations]
    unique_positions = len(set(position_keys))
    print(f"Unique position configurations: {unique_positions}")
    
    # Verify a sample
    sample_size = min(10, len(configurations))
    valid_count = 0
    for i in range(sample_size):
        if generator.verify_configuration(configurations[i]['configuration']):
            valid_count += 1
    
    print(f"Valid configurations in sample: {valid_count}/{sample_size}")


if __name__ == "__main__":
    G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
    Sampler1 = MonteCarloSampler(
        G=G, 
        led_min=1, 
        led_max=16,
        height_min=5,
        height_max=15,
        height_step=1,
        seed=42,
        sampling_mode="uniform",
        verbose=True
    )

    num_samples = 500000
    samples = Sampler1.generate_unique_samples(num_samples=num_samples, use_parallel=True)
    save_configurations(samples, Sampler1, f"results/sampled_configs/MonteCarlo_{Sampler1.led_min}-{Sampler1.led_max}_leds_{Sampler1.height_min}-{Sampler1.height_max}_cm_{num_samples}_samples.npz")

    # from geometric_model import GeometricModel
    # samples, z = GeometricModel.read_configurations(GeometricModel, "results/led_1-16_height_5-15/configurations.npz")

    # print(samples[:5])
    # print(z[:5])