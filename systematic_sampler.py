import numpy as np
from itertools import combinations
import concurrent.futures
import os
from typing import List, Tuple
from utils.grid import Grid
import json
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

class Sampler:
    """
    Generates symmetric LED configurations on a square N x N grid, placing a specified number of LEDs (M) such that each configuration satisfies both point symmetry (about the grid center) and quadrant balance constraints.

    This class is designed for systematic sampling of LED positions for irradiation modules, ensuring that the resulting configurations are physically balanced and symmetric. It supports parallelized generation of all valid configurations, visualization, and export to standard formats.

    Parameters
    ----------
    G : Grid
        Grid object defining the grid parameters (dimensions, step size, margins, etc.).
    led_count : int, optional
        Total number of LEDs to place on the grid (default: 8).
    height : int, optional
        Height of LEDs above the grid in centimeters (default: 13).
    verbose : bool, optional
        If True, prints progress and diagnostic information during configuration generation (default: True).

    Attributes
    ----------
    G : Grid
        The grid object used for configuration.
    grid_size : int
        Number of grid points along one axis (assumes square grid).
    grid_step : float
        Step size between grid points (from Grid object).
    led_count : int
        Number of LEDs to place.
    height : int
        Height of LEDs above the grid (cm).
    verbose : bool
        Controls verbosity of output.
    half_size : int
        Half the grid size (used for coordinate calculations).
    coords : list of int
        List of coordinate values for grid points.
    center_used : bool
        Whether the grid center is used (True if led_count is odd).
    pair_count : int
        Number of symmetric LED pairs (led_count // 2 or (led_count-1)//2 if odd).
    Q1_points, Q2_points, Q3_points, Q4_points : list of tuple
        Lists of grid points in each quadrant.
    pos_x_axis, pos_y_axis, all_pos_axis : list of tuple
        Lists of grid points on positive axes.
    neg_x_axis, neg_y_axis : list of tuple
        Lists of grid points on negative axes.
    center : tuple or None
        Coordinates of the grid center, or None if not applicable.

    Raises
    ------
    ValueError
        If the grid size is not odd (required for a clear center point).

    Notes
    -----
    - The grid must be square and have an odd number of points along each axis to ensure a unique center.
    - Configurations are generated such that each LED has a symmetric counterpart with respect to the grid center, except possibly one at the center if led_count is odd.
    - Quadrant balance ensures that the sum of LED contributions to each quadrant is equal.
    - The class supports parallel processing for efficient configuration generation on large grids.

    """

    def __init__(self, G: Grid, led_count: int = 8, height: int = 13, verbose: bool = True):
        """
        Initialize the symmetric LED configuration generator.

        Parameters
        ----------
        G : Grid
            Grid object defining the grid parameters (dimensions, step size, margins, etc.).
        led_count : int, optional
            Total number of LEDs to place on the grid (default: 8).
        height : int, optional
            Height of LEDs above the grid in centimeters (default: 13).
        verbose : bool, optional
            If True, prints progress and diagnostic information during configuration generation (default: True).

        Raises
        ------
        ValueError
            If the grid size is not odd (required for a clear center point).
        """

        if np.sqrt(len(G.grid_points)).astype(int) % 2 == 0:
            raise ValueError("Grid size must be odd for clear center point")

        self.G = G
        self.grid_size = np.sqrt(len(G.grid_points)).astype(int)  # Assuming square grid
        self.grid_step = G.step  
        self.led_count = led_count
        self.height = height
        self.verbose = verbose

        self.half_size = self.grid_size // 2
        self.coords = list(range(-self.half_size, self.half_size + 1))
        
        # Determine if center is used
        self.center_used = (led_count % 2 == 1)
        self.pair_count = led_count // 2  # Number of symmetric pairs
        if self.center_used:
            self.pair_count = (led_count - 1) // 2
        
        # Precompute all point classifications
        self.Q1_points, self.Q2_points, self.Q3_points, self.Q4_points = [], [], [], []
        self.pos_x_axis, self.pos_y_axis, self.all_pos_axis = [], [], []
        self.neg_x_axis, self.neg_y_axis = [], []
        self.center = (0, 0) if self.grid_size % 2 == 1 else None
        
        self._classify_points()
        
    def _classify_points(self):
        """
        Categorizes all grid points into quadrants (Q1–Q4) and axes (positive/negative x and y),
        based on their integer coordinates relative to the grid center.

        This method populates the following instance attributes:
        - ``Q1_points``: List of points in the first quadrant (x > 0, y > 0)
        - ``Q2_points``: List of points in the second quadrant (x < 0, y > 0)
        - ``Q3_points``: List of points in the third quadrant (x < 0, y < 0)
        - ``Q4_points``: List of points in the fourth quadrant (x > 0, y < 0)
        - ``pos_x_axis``: Points on the positive x-axis (x > 0, y == 0)
        - ``neg_x_axis``: Points on the negative x-axis (x < 0, y == 0)
        - ``pos_y_axis``: Points on the positive y-axis (x == 0, y > 0)
        - ``neg_y_axis``: Points on the negative y-axis (x == 0, y < 0)
        - ``all_pos_axis``: Concatenation of ``pos_x_axis`` and ``pos_y_axis``

        :raises: None
        :side effects:
            Modifies the instance attributes listed above in-place.
        :notes:
            - The grid is assumed to be square and centered at (0, 0).
            - The center point (0, 0) is excluded from all quadrant and axis lists.
            - This method should be called during initialization and whenever the grid changes.
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
        
    def get_quadrant_weight(self, point: Tuple[int, int]) -> Tuple[float, float, float, float]:
        """
        Determines the fractional contribution of a grid point to each of the four quadrants (Q1, Q2, Q3, Q4).

        The contribution is based on the point's position relative to the grid center:
        - Points in a quadrant contribute fully (1.0) to that quadrant.
        - Points on axes contribute fractionally to adjacent quadrants.
        - The center point (0, 0) contributes equally (0.25) to all quadrants.

        :param point: Tuple[int, int]
            The (x, y) coordinates of the grid point (in grid index units, not cm).
        :return: Tuple[float, float, float, float]
            The weights for quadrants Q1, Q2, Q3, Q4, respectively. Each value is in [0, 1].
        :notes:
            - The sum of the returned weights is always 1.0 for the center, and 1.0 for all other points except those not on the grid (which return all zeros).
            - Quadrant definitions:
                Q1: x > 0, y > 0
                Q2: x < 0, y > 0
                Q3: x < 0, y < 0
                Q4: x > 0, y < 0
            - Axis points split their contribution between adjacent quadrants.
        """
        x, y = point
        if x == 0 and y == 0:
            return (0.25, 0.25, 0.25, 0.25)
        elif x > 0 and y > 0:
            return (1.0, 0.0, 0.0, 0.0)
        elif x < 0 and y > 0:
            return (0.0, 1.0, 0.0, 0.0)
        elif x < 0 and y < 0:
            return (0.0, 0.0, 1.0, 0.0)
        elif x > 0 and y < 0:
            return (0.0, 0.0, 0.0, 1.0)
        elif y == 0 and x != 0:  # x-axis
            return (0.5, 0.0, 0.0, 0.5)
        elif x == 0 and y != 0:  # y-axis
            return (0.5, 0.5, 0.0, 0.0)
        return (0.0, 0.0, 0.0, 0.0)
    
    def calculate_quadrant_weights(self, configuration: List[Tuple[int, int]]) -> Tuple[float, float, float, float]:
        """
        Computes the total fractional contribution of a given LED configuration to each of the four quadrants (Q1, Q2, Q3, Q4) of the grid.

        Each LED position contributes to one or more quadrants depending on its location (quadrant, axis, or center). The sum of all contributions equals the total number of LEDs.

        :param configuration: List of (x, y) tuples
            List of LED positions, where each tuple represents integer grid coordinates (not scaled by grid step).
        :type configuration: List[Tuple[int, int]]

        :returns: Tuple containing the total weights for quadrants Q1, Q2, Q3, and Q4, respectively. Each value is a float representing the sum of contributions from all LEDs to that quadrant.
        :rtype: Tuple[float, float, float, float]

        .. note::
            The sum of the returned weights should equal the total number of LEDs in the configuration.
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
        Checks whether a given LED configuration satisfies both point symmetry (about the grid center) and quadrant balance constraints.

        The configuration is valid if:
        - For every LED at (x, y), there is also an LED at (-x, -y) (point symmetry).
        - The sum of fractional contributions to each quadrant is equal (quadrant balance).

        :param configuration: List of (x, y) tuples
            List of LED positions, where each tuple represents integer grid coordinates (not scaled by grid step).
        :type configuration: List[Tuple[int, int]]

        :returns: True if the configuration is symmetric and balanced; False otherwise.
        :rtype: bool

        .. note::
            The function assumes the configuration is on a square, odd-sized grid centered at (0, 0).

        """
        # Check point symmetry
        config_set = set(configuration)
        for point in configuration:
            symmetric_point = (-point[0], -point[1])
            if symmetric_point not in config_set:
                return False
        # Check quadrant balance
        q1, q2, q3, q4 = self.calculate_quadrant_weights(configuration)
        target_weight = self.led_count / 4  # Each quadrant should have equal weight
        return (abs(q1 - target_weight) < 1e-6 and 
                abs(q2 - target_weight) < 1e-6 and
                abs(q3 - target_weight) < 1e-6 and
                abs(q4 - target_weight) < 1e-6)

    def generate_all_configurations(self, max_workers: int | None = None) -> np.ndarray:
        """
        Systematically generates all unique, valid LED configurations on the grid that satisfy point symmetry and quadrant balance constraints.

        The method explores all possible symmetric arrangements of the specified number of LEDs, using parallel processing to accelerate computation. Each configuration is checked for symmetry and balance, and duplicates are removed. The resulting configurations are scaled by the grid step size and returned as an array.

        :param max_workers: Maximum number of worker processes to use for parallel computation. If None, uses all available CPU cores.
        :type max_workers: int or None, optional

        :returns: Array of unique valid configurations. Each entry is a numpy array of shape (N, 2), where N is the number of LEDs, and coordinates are in centimeters (scaled by grid step).
        :rtype: np.ndarray

        :raises ValueError: If the grid size is not odd (enforced during initialization).
        :side effects: May print progress information if verbosity is enabled. Uses multiprocessing and may increase system resource usage.

        .. note::
            The returned configurations are sorted, unique, and scaled to real-world units. The method may take significant time for large grids or high LED counts.
        """
        
        if self.verbose:
            print(f"Generating configurations for {self.led_count} LEDs on {self.grid_size}x{self.grid_size} grid")
            print(f"Center LED: {'Yes' if self.center_used else 'No'}")
            print(f"Symmetric pairs: {self.pair_count}")
            print(f"Target quadrant weight: {self.led_count / 4}")
        
        all_configurations = []
        
        # Find all valid (a, b, c) combinations
        valid_combinations = self._find_valid_abc_combinations()
        if self.verbose:
            print(f"Found {len(valid_combinations)} valid (a, b, c) combinations")
        
        for a, b, c in valid_combinations:
            if self.verbose:
                print(f"\nProcessing case (a={a}, b={b}, c={c})")
            if max_workers is None:  # use all available cores if not specified
                max_workers = os.cpu_count()
                
            case_configs = self._process_case_parallel(a, b, c, max_workers)
            
            all_configurations.extend(case_configs)
            print(f"Case ({a},{b},{c}): {len(case_configs)} configurations")
        
        # Remove duplicates
        unique_configs = self._remove_duplicates(all_configurations)

        unique_config_array = np.empty(len(unique_configs), dtype=object)
        for i, config in enumerate(unique_configs):
            unique_config_array[i] = np.array(config, dtype=float)

        unique_config_array = unique_config_array * self.grid_step  # Scale by grid step size

        if self.verbose:
            print(f"\n=== FINAL RESULTS ===")
            print(f"Total unique configurations: {len(unique_configs)}")
        
        return unique_config_array
    
    def _find_valid_abc_combinations(self) -> List[Tuple[int, int, int]]:
        """
        Systematically determines all valid combinations of symmetric LED pair types (a, b, c) that satisfy the quadrant balance constraint for the current grid and LED count.

        This method explores all possible ways to distribute the available symmetric LED pairs into three categories:
        - Type A: Pairs placed in the first and third quadrants (Q1/Q3)
        - Type B: Pairs placed in the second and fourth quadrants (Q2/Q4)
        - Type C: Pairs placed symmetrically on the axes
        The sum of these pairs must match the total number of symmetric pairs required for the current configuration, and the resulting arrangement must ensure that the total fractional contribution to each quadrant is equal (quadrant balance).

        :returns: List of valid (a, b, c) tuples, where:
            - a (int): Number of Type A pairs (Q1/Q3)
            - b (int): Number of Type B pairs (Q2/Q4)
            - c (int): Number of axis pairs
        :rtype: List[Tuple[int, int, int]]

        .. note::
            - The method uses the current object's ``led_count``, ``pair_count``, and ``center_used`` attributes.
            - The grid is assumed to be square and centered at (0, 0).
            - The center LED (if present) is accounted for by adjusting the target quadrant weight.
        """
        valid_combinations = []
        target_weight = self.led_count / 4
        
        # Adjust target weight for center contribution
        center_adjustment = 0.25 if self.center_used else 0.0
        
        for a in range(self.pair_count + 1):
            for b in range(self.pair_count + 1 - a):
                c = self.pair_count - a - b
                if c < 0:
                    continue
                
                # Check quadrant balance equations
                q1_weight = a + 0.5 * c + center_adjustment
                q2_weight = b + 0.5 * c + center_adjustment
                q3_weight = a + 0.5 * c + center_adjustment
                q4_weight = b + 0.5 * c + center_adjustment
                
                if (abs(q1_weight - target_weight) < 1e-6 and 
                    abs(q2_weight - target_weight) < 1e-6 and
                    abs(q3_weight - target_weight) < 1e-6 and
                    abs(q4_weight - target_weight) < 1e-6):
                    valid_combinations.append((a, b, c))
        
        return valid_combinations

    def _process_case_parallel(self, a: int, b: int, c: int, max_workers: int | None) -> List[List[Tuple[int, int]]]:
        """
        Efficiently generates all valid LED configurations for a specific (a, b, c) case using parallel processing.

        This method delegates to the appropriate internal routine depending on the values of ``a``, ``b``, and ``c``:
        - If only axis pairs are present (a == 0 and b == 0), processes axis-only configurations.
        - If no axis pairs are present (c == 0), processes configurations with only Q1/Q3 and Q2/Q4 pairs.
        - Otherwise, processes mixed cases with all three pair types.
        Parallelization is used to accelerate the generation of large numbers of configurations.

        :param a: Number of Type A symmetric pairs (Q1/Q3)
        :type a: int
        :param b: Number of Type B symmetric pairs (Q2/Q4)
        :type b: int
        :param c: Number of symmetric pairs on axes
        :type c: int
        :param max_workers: Maximum number of worker processes for parallel computation. If None, uses all available CPU cores.
        :type max_workers: int or None, optional

        :returns: List of valid configurations for the given (a, b, c) case. Each configuration is a list of (x, y) tuples (integer grid coordinates, not scaled).
        :rtype: List[List[Tuple[int, int]]]

        :side effects:
            - Uses multiprocessing (may increase system resource usage).
            - May print progress information if verbosity is enabled.

        .. note::
            - The method assumes the grid and LED count are compatible with the requested (a, b, c) values.
            - The returned configurations are not guaranteed to be unique across different (a, b, c) cases; deduplication is performed at a higher level.
        """

        if max_workers is None:
            max_workers = os.cpu_count()

        if a == 0 and b == 0:
            return self._process_axis_only_case(c, max_workers)
        elif c == 0:
            return self._process_no_axis_case(a, b, max_workers)
        else:
            return self._process_mixed_case(a, b, c, max_workers)

    def _process_axis_only_case(self, c: int, max_workers: int | None) -> List[List[Tuple[int, int]]]:
        """
        Efficiently generates all valid LED configurations for the special case where only axis pairs (i.e., symmetric pairs on the axes) are present, using parallel processing for scalability.

        This method systematically explores all possible combinations of axis pairs (from the set of positive x- and y-axis points), optionally including the grid center if required by the configuration. Each combination is expanded to a full symmetric configuration by adding both the selected axis points and their symmetric counterparts. The process is parallelized to improve performance on large grids or high axis pair counts.

        :param c: Number of axis pairs to select. Each pair consists of a point on the positive x- or y-axis and its symmetric counterpart.
        :type c: int
        :param max_workers: Maximum number of worker processes to use for parallel computation. If None, uses all available CPU cores.
        :type max_workers: int or None

        :returns:
            List of valid configurations for the axis-only case. Each configuration is a list of (x, y) tuples (integer grid coordinates, not scaled).
        :rtype: List[List[Tuple[int, int]]]

        :side effects:
            - Uses multiprocessing (may increase system resource usage).
            - May print progress information if verbosity is enabled.

        .. note::
            - The method assumes the grid and LED count are compatible with the requested number of axis pairs.
            - The returned configurations are not guaranteed to be unique across different cases; deduplication is performed at a higher level.
        """
        if c > len(self.all_pos_axis):
            return []
            
        axis_combinations = list(combinations(self.all_pos_axis, c))
        print(f"Processing {len(axis_combinations)} axis combinations")
        
        if not axis_combinations:
            return []

        chunk_size = max(1, len(axis_combinations) // (max_workers * 4)) if max_workers else len(axis_combinations)
        chunks = [axis_combinations[i:i + chunk_size] for i in range(0, len(axis_combinations), chunk_size)]
        
        all_configs = []

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self._process_axis_chunk, chunk) for chunk in chunks]
            
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                chunk_configs = future.result()
                all_configs.extend(chunk_configs)
                completed += 1
                if len(chunks) > 1 and self.verbose:
                    print(f"Axis case progress: {completed}/{len(chunks)} chunks")
        
        return all_configs
    
    def _process_axis_chunk(self, axis_chunk: List[Tuple]) -> List[List[Tuple[int, int]]]:
        """
        Processes a subset (chunk) of axis pair combinations to generate all valid symmetric LED configurations for the axis-only case.

        Each axis pair combination in the chunk is expanded to a full configuration by adding both the selected axis points and their symmetric counterparts. If the configuration requires a center LED, it is included. Only configurations with the correct total number of LEDs are returned.

        :param axis_chunk: List of axis point combinations, where each combination is a tuple of axis points (positive x- or y-axis).
        :type axis_chunk: List[Tuple[int, int]]

        :returns:
            List of valid configurations for the given chunk. Each configuration is a sorted list of (x, y) tuples (integer grid coordinates, not scaled).
        :rtype: List[List[Tuple[int, int]]]

        .. note::
            - This method is intended for internal use with parallel processing.
            - The returned configurations are not guaranteed to be unique across all chunks; deduplication is performed at a higher level.
        """
        chunk_configs = []
        for axis_points in axis_chunk:
            config = set()
            
            # Add center if needed
            if self.center_used:
                config.add(self.center)
            
            # Add axis pairs
            for axis_point in axis_points:
                config.add(axis_point)
                config.add((-axis_point[0], -axis_point[1]))
            
            if len(config) == self.led_count:
                chunk_configs.append(sorted(config))
        
        return chunk_configs

    def _process_no_axis_case(self, a: int, b: int, max_workers: int | None) -> List[List[Tuple[int, int]]]:
        """
        Systematically generates all valid LED configurations for the case where only Type A (Q1/Q3) and Type B (Q2/Q4) symmetric pairs are present, with no axis pairs included. Utilizes parallel processing to efficiently explore all possible combinations.

        This method considers all possible selections of points from the first quadrant (Q1) and second quadrant (Q2) to form symmetric pairs, ensuring that each configuration is point symmetric about the grid center and contains the correct total number of LEDs. If the configuration requires a center LED (odd LED count), it is included in each configuration. The process is parallelized by dividing the Q1 combinations into chunks, each processed independently.

        :param a: Number of Type A symmetric pairs (placed in Q1 and their symmetric counterparts in Q3).
        :type a: int
        :param b: Number of Type B symmetric pairs (placed in Q2 and their symmetric counterparts in Q4).
        :type b: int
        :param max_workers: Maximum number of worker processes for parallel computation. If None, uses all available CPU cores.
        :type max_workers: int or None, optional

        :returns:
            List of valid configurations for the given case. Each configuration is a list of (x, y) tuples (integer grid coordinates, not scaled).
        :rtype: List[List[Tuple[int, int]]]

        :side effects:
            - Uses multiprocessing (may increase system resource usage).
            - May print progress information if verbosity is enabled.

        .. note::
            - The method assumes the grid and LED count are compatible with the requested values of ``a`` and ``b``.
            - The returned configurations are not guaranteed to be unique across all cases; deduplication is performed at a higher level.
            - Each configuration is point symmetric about the grid center and contains the correct number of LEDs.
        """

        if a > len(self.Q1_points) or b > len(self.Q2_points):
            return []
        
        q1_combinations = list(combinations(self.Q1_points, a)) if a > 0 else [()]
        q2_combinations = list(combinations(self.Q2_points, b)) if b > 0 else [()]
        
        if self.verbose:
            print(f"Processing {len(q1_combinations)} Q1 combinations × {len(q2_combinations)} Q2 combinations")
        
        if not q1_combinations or not q2_combinations:
            return []

        chunk_size = max(1, len(q1_combinations) // (max_workers * 4)) if max_workers else len(q1_combinations)
        chunks = [q1_combinations[i:i + chunk_size] for i in range(0, len(q1_combinations), chunk_size)]
        
        all_configs = []

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self._process_no_axis_chunk, chunk, q2_combinations) for chunk in chunks]
            
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                chunk_configs = future.result()
                all_configs.extend(chunk_configs)
                completed += 1
                if len(chunks) > 1 and self.verbose:
                    print(f"No-axis case progress: {completed}/{len(chunks)} chunks")
        
        return all_configs
    
    def _process_no_axis_chunk(self, q1_chunk: List[Tuple], q2_combinations: List[Tuple]) -> List[List[Tuple[int, int]]]:
        """
        Processes a subset (chunk) of Q1 point combinations, pairing each with all possible Q2 combinations, to generate all valid symmetric LED configurations for the no-axis-pair case.

        For each combination of Q1 and Q2 points, constructs a configuration by adding both the selected points and their symmetric counterparts. If the configuration requires a center LED (odd LED count), it is included. Only configurations with the correct total number of LEDs are returned. This method is intended for use with parallel processing, where each chunk is processed independently.

        :param q1_chunk: List of Q1 point combinations to process in this chunk. Each entry is a tuple of (x, y) coordinates (integer grid indices).
        :type q1_chunk: List[Tuple[int, int]]
        :param q2_combinations: List of all Q2 point combinations to pair with each Q1 selection. Each entry is a tuple of (x, y) coordinates (integer grid indices).
        :type q2_combinations: List[Tuple[int, int]]

        :returns:
            List of valid configurations for the given chunk. Each configuration is a sorted list of (x, y) tuples (integer grid coordinates, not scaled).
        :rtype: List[List[Tuple[int, int]]]

        .. note::
            - This method is intended for internal use with parallel processing in the no-axis-pair case.
            - The returned configurations are not guaranteed to be unique across all chunks; deduplication is performed at a higher level.
            - Each configuration is point symmetric about the grid center and contains the correct number of LEDs.
        """
        chunk_configs = []
        for q1_selection in q1_chunk:
            for q2_selection in q2_combinations:
                config = set()
                
                # Add center if needed
                if self.center_used:
                    config.add(self.center)
                
                # Add Type A pairs (Q1 -> Q3)
                for pt in q1_selection:
                    config.add(pt)
                    config.add((-pt[0], -pt[1]))
                
                # Add Type B pairs (Q2 -> Q4)
                for pt in q2_selection:
                    config.add(pt)
                    config.add((-pt[0], -pt[1]))
                
                if len(config) == self.led_count:
                    chunk_configs.append(sorted(config))
        
        return chunk_configs

    def _process_mixed_case(self, a: int, b: int, c: int, max_workers: int | None) -> List[List[Tuple[int, int]]]:
        """
        Efficiently generates all valid LED configurations for a specific case where the configuration consists of a mix of Type A (Q1/Q3), Type B (Q2/Q4), and axis-symmetric LED pairs, using parallel processing for scalability.

        This method systematically explores all possible combinations of symmetric LED pairs distributed as follows:
        - Type A: Pairs placed in the first and third quadrants (Q1/Q3)
        - Type B: Pairs placed in the second and fourth quadrants (Q2/Q4)
        - Axis pairs: Pairs placed symmetrically on the positive axes (x or y)
        The total number of pairs is determined by the input parameters ``a``, ``b``, and ``c``. Each configuration is constructed by selecting the specified number of pairs from each category, ensuring point symmetry about the grid center and the correct total number of LEDs. If the configuration requires a center LED (odd LED count), it is included in each configuration. The process is parallelized by dividing the Q1 combinations into chunks, each processed independently.

        :param a: Number of Type A symmetric pairs (placed in Q1 and their symmetric counterparts in Q3).
        :type a: int
        :param b: Number of Type B symmetric pairs (placed in Q2 and their symmetric counterparts in Q4).
        :type b: int
        :param c: Number of symmetric pairs on the axes (each pair consists of a point on the positive x- or y-axis and its symmetric counterpart).
        :type c: int
        :param max_workers: Maximum number of worker processes for parallel computation. If None, uses all available CPU cores.
        :type max_workers: int or None, optional

        :returns:
            List of valid configurations for the given (a, b, c) case. Each configuration is a list of (x, y) tuples (integer grid coordinates, not scaled).
        :rtype: List[List[Tuple[int, int]]]

        :side effects:
            - Uses multiprocessing (may increase system resource usage).
            - May print progress information if verbosity is enabled.

        .. note::
            - The method assumes the grid and LED count are compatible with the requested values of ``a``, ``b``, and ``c``.
            - The returned configurations are not guaranteed to be unique across all cases; deduplication is performed at a higher level.
            - Each configuration is point symmetric about the grid center and contains the correct number of LEDs.
        """
        if (a > len(self.Q1_points) or b > len(self.Q2_points) or 
            c > len(self.all_pos_axis)):
            return []
        
        q1_combinations = list(combinations(self.Q1_points, a)) if a > 0 else [()]
        q2_combinations = list(combinations(self.Q2_points, b)) if b > 0 else [()]
        axis_combinations = list(combinations(self.all_pos_axis, c)) if c > 0 else [()]
        
        if self.verbose:
            print(f"Processing {len(q1_combinations)} Q1 × {len(q2_combinations)} Q2 × {len(axis_combinations)} axis")
        
        if not q1_combinations or not q2_combinations or not axis_combinations:
            return []
        
        # Split by Q1 combinations
        chunk_size = max(1, len(q1_combinations) // (max_workers * 2)) if max_workers else len(q1_combinations)
        q1_chunks = [q1_combinations[i:i + chunk_size] for i in range(0, len(q1_combinations), chunk_size)]
        
        all_configs = []

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self._process_mixed_chunk, chunk, q2_combinations, axis_combinations) 
                      for chunk in q1_chunks]
            
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                chunk_configs = future.result()
                all_configs.extend(chunk_configs)
                completed += 1
                if len(q1_chunks) > 1 and self.verbose:
                    print(f"Mixed case progress: {completed}/{len(q1_chunks)} chunks")
        
        return all_configs
    
    def _process_mixed_chunk(self, q1_chunk: List[Tuple], q2_combinations: List[Tuple], 
                           axis_combinations: List[Tuple]) -> List[List[Tuple[int, int]]]:
        """
        Generates all valid symmetric LED configurations for a subset (chunk) of possible combinations in the mixed case, where the configuration consists of a mix of Type A (Q1/Q3), Type B (Q2/Q4), and axis-symmetric LED pairs.

        This method is intended for internal use with parallel processing. For each combination of points from the provided Q1, Q2, and axis chunks, it constructs a configuration by adding both the selected points and their symmetric counterparts with respect to the grid center. If the configuration requires a center LED (odd LED count), it is included. Only configurations with the correct total number of LEDs are returned. The resulting configurations are sorted and returned as a list.

        :param q1_chunk: List of Q1 point combinations to process in this chunk. Each entry is a tuple of (x, y) coordinates (integer grid indices) representing points in the first quadrant.
        :type q1_chunk: List[Tuple[int, int]]
        :param q2_combinations: List of all Q2 point combinations to pair with each Q1 selection. Each entry is a tuple of (x, y) coordinates (integer grid indices) representing points in the second quadrant.
        :type q2_combinations: List[Tuple[int, int]]
        :param axis_combinations: List of all axis point combinations to pair with each Q1 and Q2 selection. Each entry is a tuple of (x, y) coordinates (integer grid indices) representing points on the positive axes.
        :type axis_combinations: List[Tuple[int, int]]

        :returns:
            List of valid configurations for the given chunk. Each configuration is a sorted list of (x, y) tuples (integer grid coordinates, not scaled).
        :rtype: List[List[Tuple[int, int]]]

        .. note::
            - This method is intended for internal use with parallel processing in the mixed-pair case.
            - The returned configurations are not guaranteed to be unique across all chunks; deduplication is performed at a higher level.
            - Each configuration is point symmetric about the grid center and contains the correct number of LEDs.
            - The method assumes the grid and LED count are compatible with the requested values of ``a``, ``b``, and ``c``.
        """
        chunk_configs = []
        for q1_selection in q1_chunk:
            for q2_selection in q2_combinations:
                for axis_selection in axis_combinations:
                    config = set()
                    
                    # Add center if needed
                    if self.center_used:
                        config.add(self.center)
                    
                    # Type A pairs
                    for pt in q1_selection:
                        config.add(pt)
                        config.add((-pt[0], -pt[1]))
                    
                    # Type B pairs
                    for pt in q2_selection:
                        config.add(pt)
                        config.add((-pt[0], -pt[1]))
                    
                    # Axis pairs
                    for axis_point in axis_selection:
                        config.add(axis_point)
                        config.add((-axis_point[0], -axis_point[1]))
                    
                    if len(config) == self.led_count:
                        chunk_configs.append(sorted(config))
        
        return chunk_configs
    
    def _remove_duplicates(self, configurations: List[List[Tuple[int, int]]]) -> List[List[Tuple[int, int]]]:
        """
        Removes duplicate LED configurations from a list, ensuring each configuration is unique.

        This method compares each configuration (a list of (x, y) tuples representing LED positions) and eliminates duplicates by converting each configuration to a tuple of tuples and storing it in a set for fast lookup. Only the first occurrence of each unique configuration is retained in the output list. The order of configurations in the returned list corresponds to their first appearance in the input.

        :param configurations: List of configurations to process. Each configuration is a list of (x, y) tuples, where each tuple represents the integer grid coordinates of an LED position (not scaled by grid step).
        :type configurations: List[List[Tuple[int, int]]]

        :returns: List of unique configurations, with duplicates removed. Each configuration is a list of (x, y) tuples as in the input.
        :rtype: List[List[Tuple[int, int]]]

        .. note::
            - Configurations are considered duplicates if their sorted list of (x, y) tuples matches exactly.
            - The method does not sort the input list; it preserves the order of first occurrence.
            - This method is typically used after generating all possible configurations to ensure uniqueness before further processing or saving.
        """
        seen = set()
        unique_configs = []
        
        for config in configurations:
            config_tuple = tuple(tuple(point) for point in config)
            if config_tuple not in seen:
                seen.add(config_tuple)
                unique_configs.append(config)
        
        return unique_configs
    
    def plot_configuration(self, ax: Axes, configuration: List[Tuple[int, int]], z: float = -1., point_size: int = 20, set_title: bool = False) -> Axes:
        """
        Visualizes a given LED configuration on a 2D grid using Matplotlib, displaying the LED positions and grid layout.

        This method plots the provided configuration of LED positions on the specified Matplotlib Axes object. The grid is scaled and labeled according to the underlying Grid object. Optionally, the plot can include a title indicating the number of LEDs and their height above the grid.

        :param ax: Matplotlib Axes
            The Axes object on which to plot the configuration.
        :type ax: matplotlib.axes.Axes
        :param configuration: List of (x, y) tuples
            The LED positions to plot, as a list of tuples representing integer grid coordinates (not scaled by grid step).
        :type configuration: List[Tuple[int, int]]
        :param z: Height of LEDs above the grid in centimeters, used for the plot title. If set to -1 (default), uses the Sampler's height attribute.
        :type z: float, optional
        :param point_size: Size of the plotted LED points.
        :type point_size: int, optional
        :param set_title: Whether to set a title on the plot indicating configuration details.
        :type set_title: bool, optional

        :return: The Matplotlib Axes object with the plotted configuration.
        :rtype: matplotlib.axes.Axes

        :side effects:
            - Modifies the provided Axes object by adding scatter points, grid lines, and labels.
            - Sets Matplotlib styles for consistent appearance.

        .. note::
            - The configuration coordinates are interpreted as grid indices and are not automatically scaled to centimeters; scaling should be handled prior to plotting if needed.
            - The grid limits and ticks are determined by the associated Grid object.
            - The method does not display the plot; it only modifies the Axes object.
        """

        # set style
        plt.style.use("ICIWstyle")
        plt.style.use("visualization/publication_style.mplstyle")

        x_coords, y_coords = zip(*configuration)
        # convert to float for plotting
        x_coords = [float(x) for x in x_coords]
        y_coords = [float(y) for y in y_coords]

        ax.scatter(x_coords, y_coords, c='blue', s=point_size)

        ax.set_xlim(-self.G.width/2, self.G.width/2)
        ax.set_ylim(-self.G.height/2, self.G.height/2)
        ax.set_xticks(np.arange(-self.G.width/2 + self.G.side_margin, self.G.width/2 - self.G.side_margin + self.G.step, self.G.step), minor=True)
        ax.set_yticks(np.arange(-self.G.height/2 + self.G.top_bottom_margin, self.G.height/2 - self.G.top_bottom_margin + self.G.step, self.G.step), minor=True)
        ax.tick_params(axis='both', which='both', direction='in')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.set_aspect('equal')
        if set_title:
            ax.set_title(f"Configuration with {len(configuration)} LEDs at height {z if z != -1. else self.height} cm")

        ax.set_xlabel("$X$ / cm")
        ax.set_ylabel("$Y$ / cm")

        return ax

    def plot_configuration_small(self, ax: Axes, configuration: List[Tuple[int, int]], z: float = -1., point_size: int = 10) -> Axes:
        """
        Creates a compact visualization of an LED configuration on a Matplotlib Axes, suitable for use in subplots or summary figures.

        This method plots the provided LED positions on the given Axes object, using smaller point sizes and minimal axis labeling for a cleaner, more compact appearance. The grid is drawn according to the associated Grid object, but axis labels and tick labels are hidden to save space.

        :param ax: Matplotlib Axes
            The Axes object on which to plot the configuration.
        :type ax: matplotlib.axes.Axes
        :param configuration: List of (x, y) tuples
            The LED positions to plot, as a list of tuples representing integer grid coordinates (not scaled by grid step).
        :type configuration: List[Tuple[int, int]]
        :param z: Height of LEDs above the grid in centimeters, used for the plot title (not shown in this small plot). Default is -1.
        :type z: float, optional
        :param point_size: Size of the plotted LED points (smaller than in the main plot).
        :type point_size: int, optional

        :return: The Matplotlib Axes object with the plotted configuration.
        :rtype: matplotlib.axes.Axes

        :side effects:
            - Modifies the provided Axes object by adding scatter points, grid lines, and hiding axis labels and tick labels.
            - Sets Matplotlib styles for consistent appearance.

        .. note::
            - This method is intended for use in multi-panel figures or summary plots where space is limited.
            - The configuration coordinates are interpreted as grid indices and are not automatically scaled to centimeters.
            - The plot does not display axis labels or tick labels, and the grid is shown behind the points.
        """

        plt.style.use("ICIWstyle")
        plt.style.use("visualization/publication_style.mplstyle")

        x_coords, y_coords = zip(*configuration)
        # convert to float for plotting
        x_coords = [float(x) for x in x_coords]
        y_coords = [float(y) for y in y_coords]

        ax.scatter(x_coords, y_coords, s=point_size)
        ax.set_xlim(-self.G.width/2, self.G.width/2)
        ax.set_ylim(-self.G.height/2, self.G.height/2)
        ax.set_xticks(np.arange(-self.G.width/2 + self.G.side_margin, self.G.width/2 - self.G.side_margin + self.G.step, self.G.step), minor=True)
        ax.set_yticks(np.arange(-self.G.height/2 + self.G.top_bottom_margin, self.G.height/2 - self.G.top_bottom_margin + self.G.step, self.G.step), minor=True)
        ax.tick_params(axis='both', which='both', direction='in')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.set_aspect('equal')

        # also show ticks on top and right
        ax.tick_params(axis='x', which='both', top=True)
        ax.tick_params(axis='y', which='both', right=True)
        # dont show ticks and labels
        ax.set_xticklabels([], minor=False)
        ax.set_yticklabels([], minor=False)
        ax.set_xlabel("")
        ax.set_ylabel("")

        # hide grid behind the points
        ax.set_axisbelow(True)
        
        return ax
    
    def save_configuration_to_json(self, configuration: List[Tuple[int, int]], height: int, filename: str, save_fig_type: str = "") -> None:
        """
        Serializes a single LED configuration to a JSON file, optionally saving a visualization of the configuration as an image file.

        The method writes the configuration, including the number of LEDs, their (x, y) positions (as provided), and the LED height, to a JSON file at the specified path. If a file type is specified via ``save_fig_type``, a plot of the configuration is also saved in the same location with the corresponding extension.

        :param configuration: List of (x, y) tuples
            The LED positions to save, as a list of tuples representing grid coordinates (not scaled by grid step).
        :type configuration: List[Tuple[int, int]]
        :param height: Height of LEDs above the grid in centimeters.
        :type height: int
        :param filename: Path to the output JSON file.
        :type filename: str
        :param save_fig_type: File type for saving a plot of the configuration (e.g., 'png', 'svg'). If empty (default), no plot is saved.
        :type save_fig_type: str, optional

        :side effects:
            - Writes a JSON file to disk containing the configuration data.
            - If ``save_fig_type`` is provided, creates and saves a plot image file in the same directory.
            - May create or overwrite files at the specified paths.

        .. note::
            - The configuration is saved as provided; coordinates are not scaled to centimeters unless done prior to calling this method.
            - The plot image, if saved, uses the ``plot_configuration`` method for visualization.
        """

        with open(filename, 'w') as f:
            json.dump({
                "Number of LEDs": len(configuration),
                "(x, y) Positions in [cm]": configuration,
                "Height in [cm]": int(height),
            }, f, indent=2)
            
        if save_fig_type != "":
            fig, ax = plt.subplots(figsize=(4,4))
            ax = self.plot_configuration(ax, configuration, height, set_title=True)
            fig.savefig(filename.replace('.json', f'.{save_fig_type}'), dpi=300)
            plt.close(fig)
            
    def read_configuration_from_json(self, filename: str) -> Tuple[List[Tuple[int, int]], int]:
        """
        Loads a single LED configuration and its associated height from a JSON file.

        The JSON file is expected to contain the following keys:
        - "(x, y) Positions in [cm]": List of LED positions as (x, y) coordinate pairs (in centimeters).
        - "Height in [cm]": Integer or float specifying the height of the LEDs above the grid (in centimeters).

        :param filename: Path to the input JSON file containing the configuration.
        :type filename: str

        :returns:
            - configuration (List[Tuple[int, int]]): List of (x, y) tuples representing LED positions in centimeters.
            - height (int): Height of LEDs above the grid in centimeters.

        :raises FileNotFoundError: If the specified file does not exist.
        :raises KeyError: If required keys are missing from the JSON file.
        :raises json.JSONDecodeError: If the file is not valid JSON.

        :side effects:
            - Opens and reads from disk.

        .. note::
            The coordinates are returned as tuples, and the height is returned as an integer or float depending on the file content.
        """
        import json
        with open(filename, 'r') as f:
            data = json.load(f)
        configuration = [tuple(pos) for pos in data["(x, y) Positions in [cm]"]]
        height = data["Height in [cm]"]
        return configuration, height


def save_configurations(configurations: np.ndarray,
                        generator: Sampler, 
                        filename: str, 
                        ):
    """
    Saves a collection of LED configurations and associated metadata to a compressed NumPy NPZ file.

    The output file contains the following arrays and metadata:
        - 'configurations': Array of configurations (each configuration is an array of LED positions, typically shape (N, 2)).
        - 'total_count': Total number of configurations saved (int).
        - 'grid_size': Number of grid points along one axis (int).
        - 'led_count': Number of LEDs per configuration (int).
        - 'height': Height of LEDs above the grid in centimeters (int or float).
        - 'generator': String identifier for the generator used.
        - 'description': Human-readable description of the dataset.

    :param configurations: Array of configurations to save. Each entry should be an array-like of shape (led_count, 2), with coordinates in centimeters.
    :type configurations: np.ndarray
    :param generator: Sampler object used to generate the configurations. Used for metadata extraction.
    :type generator: Sampler
    :param filename: Path to the output NPZ file.
    :type filename: str

    :returns: None

    :raises IOError: If the file cannot be written.

    :side effects:
        - Writes a compressed NPZ file to disk at the specified location.
        - Prints a confirmation message to stdout upon successful save.

    .. note::
        The configurations are saved as provided; ensure they are scaled to the desired units (e.g., centimeters) before calling this function.
    """
    np.savez_compressed(
        filename,
        configurations=configurations,
        total_count=len(configurations),
        grid_size=generator.grid_size,
        led_count=generator.led_count,
        height=generator.height,
        generator="SymmetricLEDGenerator",
        description=f"{len(configurations)} configurations of {generator.led_count} LEDs on {generator.grid_size}x{generator.grid_size} grid (stepsize {generator.grid_step}) with height {generator.height}. Generated by SymmetricLEDGenerator."
    )
    print(f"Saved {len(configurations)} configurations to {filename}")


if __name__ == "__main__":

    G = Grid(width=33, height=34, step=2.5, side_space=1.5, top_bottom_space=2.)
    Sampler1 = Sampler(G=G, led_count=8, height=13)
    
    # generate all configurations
    configurations = Sampler1.generate_all_configurations()
    save_configurations(configurations, Sampler1, f"results/sampled_configs/systematic_{Sampler1.led_count}_leds_{Sampler1.height}_cm_{len(configurations)}_samples.npz")
    
    # plot and save best manual configuration as example
    # conf = [
    #     (-7.5,0), (7.5,0), (-10, -10), (10, 10), (-10, 10), (10, -10), (0,12.5), (0,-12.5)
    # ]
    # fig, ax = plt.subplots(figsize=(2,2))
    # ax = Sampler1.plot_configuration(ax, conf, z=13, point_size=20)
    # plt.show()

    # Sampler1.save_configuration_to_json(
    #     conf,
    #     height=13,
    #     filename="results/sampled_configs/best_manual.json",
    #     save_fig_type="svg"
    # )