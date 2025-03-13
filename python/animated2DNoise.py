import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.animation as animation
from matplotlib.widgets import Button, RadioButtons, Slider

# Import C++ bindings
import noise_demo as nd

class TransferFunctionVisualizer:
    def __init__(self):
        # Create terrain generator
        self.terrain_gen = nd.TerrainGen2DAnimation(seed=1354)
        self.current_frame = 0
        self.width, self.height = 1000, 1000  # Reduced size for faster generation
        
        # Initialize parameters
        self.params = self.create_terrain_params()
        
        # Generate the terrain data
        print("Generating animated terrain...")
        self.terrain = self.terrain_gen.generateAnimatedTerrain(self.width, self.height, self.params)
        self.frame_names = self.terrain_gen.getAnimationFrameNames(self.params)
        
        # Extract data arrays
        print("Extracting data arrays...")
        self.extract_terrain_data()
        
        # Setup visualization
        self.setup_visualization()

    def create_terrain_params(self):
        # Define parameters with transfer functions
        params = nd.TerrainParams()
        
        # Base noise parameters
        params.height_octaves = 6
        params.height_persistence = 0.5
        params.height_lacunarity = 2.0
        params.height_scale = 0.005
        
        # Continentalness parameters
        params.continental_scale = 0.001
        params.continental_octaves = 3
        
        # Ridge formation parameters
        params.ridge_scale = 0.003
        params.ridge_octaves = 2
        params.ridge_persistence = 0.5
        params.ridge_influence = 0.6
        
        # Erosion parameters
        params.erosion_scale = 0.01
        params.erosion_octaves = 4
        
        # Temperature parameters
        params.temperature_scale = 0.008
        params.temp_latitude_influence = 0.3
        
        # Humidity parameters
        params.humidity_scale = 0.004
        params.humidity_octaves = 3
        
        # Animation settings
        params.animation_frames = 8
        
        # Turn on auto-generation of transfer functions
        params.auto_generate_transfer_functions = True
        
        # Disable old spline auto-generation
        params.auto_generate_splines = False
        
        """
        # EXAMPLE: Custom transfer functions for each parameter
        # Uncomment and modify these to customize terrain generation
        
        # 1. Continentalness Transfer Function
        # Creates distinct continental shelves and ocean depth zones
        continental_transfer = nd.create_transfer_spline(
            "Continental Shelf Transfer", 
            nd.TerrainParameter.Continentalness,
            tension=0.5,
            weight=1.0
        )
        
        # Add control points to define the transfer curve (input → output mapping)
        # Deep ocean values
        nd.add_transfer_control_point(continental_transfer, -1.0, -1.0)   # Deepest ocean
        nd.add_transfer_control_point(continental_transfer, -0.8, -0.85)  # Deep ocean
        
        # Continental shelf - notice the plateau effect
        nd.add_transfer_control_point(continental_transfer, -0.7, -0.35)  # Shelf start
        nd.add_transfer_control_point(continental_transfer, -0.4, -0.3)   # Shelf plateau
        
        # Transition to land
        nd.add_transfer_control_point(continental_transfer, -0.2, -0.1)   # Near shore
        nd.add_transfer_control_point(continental_transfer, 0.0, 0.0)     # Coastline
        
        # Land elevation gradient
        nd.add_transfer_control_point(continental_transfer, 0.3, 0.3)     # Lowlands
        nd.add_transfer_control_point(continental_transfer, 0.6, 0.7)     # Highlands
        nd.add_transfer_control_point(continental_transfer, 1.0, 1.0)     # Maximum height
        
        # 2. Ridge Formation Transfer Function
        # Enhances mountain ridges while suppressing low values
        ridge_transfer = nd.create_transfer_spline(
            "Ridge Enhancement",
            nd.TerrainParameter.RidgeFormation,
            tension=0.3,
            weight=1.0
        )
        
        # Define the ridge transfer curve - make ridges more dramatic
        nd.add_transfer_control_point(ridge_transfer, -1.0, -0.3)     # Minimize negative values
        nd.add_transfer_control_point(ridge_transfer, -0.5, -0.1)
        nd.add_transfer_control_point(ridge_transfer, 0.0, 0.0)       # Neutral point
        nd.add_transfer_control_point(ridge_transfer, 0.3, 0.5)       # Start enhancing
        nd.add_transfer_control_point(ridge_transfer, 0.6, 0.8)       # Enhance mid-high ridges
        nd.add_transfer_control_point(ridge_transfer, 0.8, 0.95)      # Strong peaks
        nd.add_transfer_control_point(ridge_transfer, 1.0, 1.0)       # Maximum
        
        # 3. Erosion Transfer Function
        # Create more nuanced erosion effects
        erosion_transfer = nd.create_transfer_spline(
            "Erosion Pattern",
            nd.TerrainParameter.Erosion,
            tension=0.5,
            weight=1.0
        )
        
        # S-curve gives more pronounced erosion in mid-ranges
        nd.add_transfer_control_point(erosion_transfer, -1.0, -0.8)
        nd.add_transfer_control_point(erosion_transfer, -0.7, -0.6)
        nd.add_transfer_control_point(erosion_transfer, -0.3, -0.1)
        nd.add_transfer_control_point(erosion_transfer, 0.0, 0.0)
        nd.add_transfer_control_point(erosion_transfer, 0.3, 0.1)
        nd.add_transfer_control_point(erosion_transfer, 0.7, 0.6)
        nd.add_transfer_control_point(erosion_transfer, 1.0, 0.8)
        
        # 4. Temperature Transfer Function
        # Create distinct temperature bands for clearer climate zones
        temp_transfer = nd.create_transfer_spline(
            "Temperature Bands",
            nd.TerrainParameter.Temperature,
            tension=0.3,
            weight=1.0
        )
        
        # Define temperature bands with plateaus
        nd.add_transfer_control_point(temp_transfer, -1.0, -1.0)      # Polar
        nd.add_transfer_control_point(temp_transfer, -0.8, -0.8)
        nd.add_transfer_control_point(temp_transfer, -0.6, -0.6)      # Subpolar
        nd.add_transfer_control_point(temp_transfer, -0.4, -0.3)
        nd.add_transfer_control_point(temp_transfer, -0.2, -0.1)      # Temperate
        nd.add_transfer_control_point(temp_transfer, 0.0, 0.1)
        nd.add_transfer_control_point(temp_transfer, 0.2, 0.3)        # Subtropical
        nd.add_transfer_control_point(temp_transfer, 0.5, 0.6)
        nd.add_transfer_control_point(temp_transfer, 0.7, 0.8)        # Tropical
        nd.add_transfer_control_point(temp_transfer, 1.0, 1.0)
        
        # 5. Humidity Transfer Function
        # Create more diverse moisture patterns
        humidity_transfer = nd.create_transfer_spline(
            "Humidity Distribution",
            nd.TerrainParameter.Humidity,
            tension=0.4,
            weight=1.0
        )
        
        # Define humidity zones with plateaus for desert, moderate, and rainforest regions
        nd.add_transfer_control_point(humidity_transfer, -1.0, -1.0)   # Extremely dry
        nd.add_transfer_control_point(humidity_transfer, -0.7, -0.9)   # Desert plateau
        nd.add_transfer_control_point(humidity_transfer, -0.4, -0.7)
        nd.add_transfer_control_point(humidity_transfer, -0.2, -0.2)   # Semi-arid transition
        nd.add_transfer_control_point(humidity_transfer, 0.0, 0.0)     # Moderate
        nd.add_transfer_control_point(humidity_transfer, 0.2, 0.2)     # Moderate-humid
        nd.add_transfer_control_point(humidity_transfer, 0.5, 0.7)     # Humid
        nd.add_transfer_control_point(humidity_transfer, 0.8, 0.9)     # Very humid plateau
        nd.add_transfer_control_point(humidity_transfer, 1.0, 1.0)     # Maximum humidity
        
        # Add the transfer functions to the parameters
        nd.add_transfer_function(params, nd.TerrainParameter.Continentalness, continental_transfer)
        nd.add_transfer_function(params, nd.TerrainParameter.RidgeFormation, ridge_transfer)
        nd.add_transfer_function(params, nd.TerrainParameter.Erosion, erosion_transfer)
        nd.add_transfer_function(params, nd.TerrainParameter.Temperature, temp_transfer)
        nd.add_transfer_function(params, nd.TerrainParameter.Humidity, humidity_transfer)
        """
        
        return params

    def extract_terrain_data(self):
        # Extract data from terrain cells into NumPy arrays
        self.height_data = nd.extract_height_data(self.terrain)
        self.base_noise_data = nd.extract_base_noise_data(self.terrain)
        
        # Extract both raw and transformed parameter values
        # Raw values (before transfer function)
        self.raw_temperature_data = nd.extract_raw_temperature_data(self.terrain)
        self.raw_humidity_data = nd.extract_raw_humidity_data(self.terrain)
        self.raw_erosion_data = nd.extract_raw_erosion_data(self.terrain)
        self.raw_continentalness_data = nd.extract_raw_continentalness_data(self.terrain)
        self.raw_ridge_data = nd.extract_raw_ridge_value_data(self.terrain)
        
        # Transformed values (after transfer function)
        self.temperature_data = nd.extract_temperature_data(self.terrain)
        self.humidity_data = nd.extract_humidity_data(self.terrain)
        self.erosion_data = nd.extract_erosion_data(self.terrain)
        self.continentalness_data = nd.extract_continentalness_data(self.terrain)
        self.ridge_data = nd.extract_ridge_value_data(self.terrain)
        
        self.biome_data = nd.extract_biome_data(self.terrain)
        
        # Scale data for visualization - convert from [-1,1] to [0,1] range
        self.height_scaled = (self.height_data + 1) / 2
        self.base_noise_scaled = (self.base_noise_data + 1) / 2
        
        # Raw data scaling
        self.raw_temperature_scaled = (self.raw_temperature_data + 1) / 2
        self.raw_humidity_scaled = (self.raw_humidity_data + 1) / 2
        self.raw_erosion_scaled = (self.raw_erosion_data + 1) / 2
        self.raw_continentalness_scaled = (self.raw_continentalness_data + 1) / 2
        self.raw_ridge_scaled = (self.raw_ridge_data + 1) / 2
        
        # Transformed data scaling
        self.temperature_scaled = (self.temperature_data + 1) / 2
        self.humidity_scaled = (self.humidity_data + 1) / 2
        self.erosion_scaled = (self.erosion_data + 1) / 2
        self.continentalness_scaled = (self.continentalness_data + 1) / 2
        self.ridge_scaled = (self.ridge_data + 1) / 2
        
        # Create biome color map with expanded biome list
        self.biome_colors = {
            # Water biomes
            int(nd.BiomeType.DeepOcean): (0.0, 0.0, 0.6),         # Deep blue
            int(nd.BiomeType.Ocean): (0.0, 0.0, 0.8),             # Medium blue
            int(nd.BiomeType.WarmOcean): (0.0, 0.2, 0.9),         # Lighter blue with hint of turquoise
            int(nd.BiomeType.FrozenOcean): (0.7, 0.8, 0.9),       # Light icy blue
            
            # Coastal biomes
            int(nd.BiomeType.Beach): (0.95, 0.95, 0.6),           # Light sand color
            int(nd.BiomeType.RockyShore): (0.6, 0.6, 0.6),        # Gray
            int(nd.BiomeType.Marsh): (0.4, 0.5, 0.3),             # Murky green-brown
            
            # Hot dry biomes
            int(nd.BiomeType.Desert): (0.95, 0.8, 0.5),           # Light tan
            int(nd.BiomeType.Savanna): (0.85, 0.75, 0.4),         # Darker tan
            int(nd.BiomeType.Shrubland): (0.7, 0.65, 0.3),        # Brownish-green
            int(nd.BiomeType.TropicalShrubland): (0.65, 0.55, 0.2), # Deeper brownish-green
            
            # Grassland biomes
            int(nd.BiomeType.Plains): (0.6, 0.8, 0.4),            # Light green
            int(nd.BiomeType.TemperateGrassland): (0.5, 0.75, 0.3), # Medium green
            int(nd.BiomeType.BorealPlains): (0.5, 0.7, 0.5),      # Blueish-green
            
            # Cold flat biomes
            int(nd.BiomeType.Tundra): (0.7, 0.7, 0.7),            # Light gray
            
            # Tropical forest biomes
            int(nd.BiomeType.TropicalRainforest): (0.0, 0.5, 0.0), # Deep green
            int(nd.BiomeType.TropicalSeasonalForest): (0.2, 0.6, 0.2), # Slightly lighter green
            
            # Temperate forest biomes
            int(nd.BiomeType.TemperateRainforest): (0.1, 0.6, 0.3), # Green with blue tint
            int(nd.BiomeType.DeciduousForest): (0.3, 0.65, 0.25),  # Bright green
            int(nd.BiomeType.Forest): (0.2, 0.6, 0.1),            # Standard forest green
            int(nd.BiomeType.MixedForest): (0.3, 0.55, 0.2),      # Varied green
            
            # Boreal forest biomes
            int(nd.BiomeType.Taiga): (0.2, 0.5, 0.4),             # Dark green with blue tint
            int(nd.BiomeType.SnowyTaiga): (0.4, 0.6, 0.5),        # Lighter green-blue
            
            # Elevated terrain
            int(nd.BiomeType.Hills): (0.5, 0.5, 0.4),             # Light brown
            int(nd.BiomeType.Badlands): (0.7, 0.4, 0.2),          # Reddish-brown
            int(nd.BiomeType.Mountain): (0.5, 0.5, 0.5),          # Gray
            int(nd.BiomeType.SnowyMountain): (0.7, 0.7, 0.8),     # Light blue-gray
            int(nd.BiomeType.SnowyCap): (1.0, 1.0, 1.0),          # White
            
            # Ice biomes
            int(nd.BiomeType.IceSheet): (0.8, 0.9, 1.0),          # Pale blue
        }
        
        self.biome_cmap_list = [self.biome_colors[i] for i in range(len(self.biome_colors))]
        self.biome_cmap = LinearSegmentedColormap.from_list("biome_cmap", self.biome_cmap_list, N=len(self.biome_colors))

    def setup_visualization(self):
        print("Setting up visualization...")
        # Create figure with specific layout
        self.fig = plt.figure(figsize=(16, 12))
        
        # Use GridSpec for a more flexible layout
        # Main grid divides the figure into 4 rows
        gs_main = plt.GridSpec(4, 2, height_ratios=[3, 3, 2, 1], width_ratios=[3, 1], figure=self.fig)
        
        # Create main plot for terrain animation (left side)
        self.main_ax = self.fig.add_subplot(gs_main[0:2, 0])
        self.main_ax.set_title("Terrain Development", fontsize=16)
        self.main_ax.axis('off')
        
        # Create biome plot (right side, same height as main plot)
        self.biome_ax = self.fig.add_subplot(gs_main[0:2, 1])
        self.biome_ax.set_title("Final Biomes", fontsize=16)
        self.biome_ax.axis('off')
        
        # Create parameter comparison grid (shows raw vs transformed values)
        gs_params = plt.GridSpec(2, 3, bottom=0.3, top=0.45, left=0.05, right=0.95, figure=self.fig)
        
        # Transfer function visualization
        gs_transfer = plt.GridSpec(1, 1, bottom=0.05, top=0.25, left=0.1, right=0.9, figure=self.fig)
        self.transfer_ax = self.fig.add_subplot(gs_transfer[0, 0])
        self.transfer_ax.set_title("Transfer Function", fontsize=12)
        self.transfer_ax.set_xlabel("Input Value")
        self.transfer_ax.set_ylabel("Output Value")
        self.transfer_ax.grid(True)
        
        # Parameter names for visualization
        self.parameter_names = ["Continentalness", "Ridge Formation", "Erosion", 
                          "Temperature", "Humidity"]
        
        # Create axes for raw and transformed data comparison
        self.param_axes = []
        param_titles = ["Raw vs. Transformed: " + name for name in self.parameter_names]
        
        # Create 2x3 grid of parameter visualizations (raw vs transformed)
        for i in range(5):
            row = i // 3
            col = i % 3
            ax = self.fig.add_subplot(gs_params[row, col])
            ax.set_title(param_titles[i], fontsize=10)
            ax.axis('off')
            self.param_axes.append(ax)
        
        # Create colormaps
        self.height_cmap = plt.cm.terrain
        self.temp_cmap = plt.cm.coolwarm
        self.humidity_cmap = plt.cm.Blues
        self.erosion_cmap = plt.cm.Reds
        self.continental_cmap = plt.cm.YlGn
        self.ridge_cmap = plt.cm.plasma
        
        # Initialize main image with first frame
        self.main_img = self.main_ax.imshow(
            self.height_scaled[:, :, 0], 
            cmap=self.height_cmap, 
            vmin=0, vmax=1
        )
        
        # Initialize biome image (always shows final biome state)
        self.biome_img = self.biome_ax.imshow(
            self.biome_data[:, :, -1] if self.biome_data.shape[2] > 0 else np.zeros_like(self.height_data[:, :, 0]), 
            cmap=self.biome_cmap, 
            vmin=0, vmax=len(self.biome_colors)-1
        )
        
        # Create colorbar for biomes
        cb_biome = plt.colorbar(
            self.biome_img, 
            ax=self.biome_ax, 
            fraction=0.046, 
            pad=0.04, 
            ticks=range(len(self.biome_colors))
        )
        cb_biome.set_ticklabels([
            'DeepOcean', 'Ocean', 'WarmOcean', 'FrozenOcean',
            'Beach', 'RockyShore', 'Marsh',
            'Desert', 'Savanna', 'Shrubland', 'TropicalShrubland',
            'Plains', 'TempGrassland', 'BorealPlains',
            'Tundra',
            'TropRainforest', 'TropSeasonalForest',
            'TempRainforest', 'DeciduousForest', 'Forest', 'MixedForest',
            'Taiga', 'SnowyTaiga',
            'Hills', 'Badlands', 'Mountain', 'SnowyMountain', 'SnowyCap',
            'IceSheet'
        ])
        cb_biome.ax.tick_params(labelsize=6)  # Reduce font size for the expanded list
        
        # We'll set up the parameter comparison visuals in update_frame
        self.param_imgs_raw = []
        self.param_imgs_transformed = []
        
        # Setup parameter selection
        ax_radio = plt.axes([0.02, 0.1, 0.08, 0.15])
        self.parameter_selector = RadioButtons(
            ax_radio, 
            ['Continentalness', 'Ridge', 'Erosion', 'Temperature', 'Humidity'],
            active=0
        )
        self.parameter_selector.on_clicked(self.on_parameter_select)
        
        # Setup transfer function plot with default continentalness
        self.selected_parameter = nd.TerrainParameter.Continentalness
        self.update_transfer_curve()
        
        # Create navigation controls
        self.setup_controls()
        
        # Descriptions for each frame
        self.descriptions = [
            "Step 1: Base Noise - Raw noise patterns",
            "Step 2: Continentalness - Transfer function applied to raw values",
            "Step 3: Ridge Formation - Transfer function transforms ridge values",
            "Step 4: Erosion - Transfer function adjusts erosion effects",
            "Step 5: Temperature - Climate zones with transfer functions",
            "Step 6: Humidity - Moisture patterns from transfer functions",
            "Step 7: Final Height - All transformed parameters combined",
            "Step 8: Biome Classification - Environmental zones based on parameters"
        ]
        
        # Main title
        self.main_title = self.fig.suptitle(
            self.descriptions[0], 
            fontsize=16, 
            y=0.95
        )
        
        # Initial display
        self.update_frame(0)
        
        # Adjust layout for better spacing
        plt.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.05, wspace=0.3, hspace=0.4)
        
    def setup_controls(self):
        # Add buttons for previous/next frame
        ax_prev = plt.axes([0.35, 0.02, 0.1, 0.03])
        ax_next = plt.axes([0.55, 0.02, 0.1, 0.03])
        
        self.btn_prev = Button(ax_prev, 'Previous')
        self.btn_next = Button(ax_next, 'Next')
        
        self.btn_prev.on_clicked(self.prev_frame)
        self.btn_next.on_clicked(self.next_frame)
        
    def update_transfer_curve(self):
        # Clear the transfer function plot
        self.transfer_ax.clear()
        
        # Get default transfer function for the selected parameter
        transfer_spline = nd.TerrainGen2DAnimation.generateDefaultTransferFunction(self.selected_parameter)[0]
        
        # Generate curve points
        curve_points = nd.generate_transfer_curve(transfer_spline, 100)
        
        # Extract x and y values
        x_values = [point[0] for point in curve_points]
        y_values = [point[1] for point in curve_points]
        
        # Extract control points for markers
        control_x = [point.inputValue for point in transfer_spline.controlPoints]
        control_y = [point.outputValue for point in transfer_spline.controlPoints]
        
        # Plot transfer curve
        self.transfer_ax.plot(x_values, y_values, '-', color='blue', label='Transfer Function')
        
        # Plot control points
        self.transfer_ax.plot(control_x, control_y, 'o', color='red', label='Control Points')
        
        # Plot reference line (y = x)
        self.transfer_ax.plot([-1, 1], [-1, 1], '--', color='gray', alpha=0.5, label='Identity')
        
        # Set limits and labels
        self.transfer_ax.set_xlim(-1.1, 1.1)
        self.transfer_ax.set_ylim(-1.1, 1.1)
        self.transfer_ax.set_title(f"Transfer Function: {self.parameter_names[int(self.selected_parameter) - 1]}", fontsize=12)
        self.transfer_ax.legend(loc='upper left')
        self.transfer_ax.grid(True)
        
        # Redraw
        self.fig.canvas.draw_idle()
        
    def on_parameter_select(self, label):
        # Map the label to the TerrainParameter enum
        param_map = {
            'Continentalness': nd.TerrainParameter.Continentalness,
            'Ridge': nd.TerrainParameter.RidgeFormation,
            'Erosion': nd.TerrainParameter.Erosion,
            'Temperature': nd.TerrainParameter.Temperature,
            'Humidity': nd.TerrainParameter.Humidity
        }
        
        self.selected_parameter = param_map[label]
        self.update_transfer_curve()
        
    def update_frame(self, frame_idx):
        self.current_frame = frame_idx
        
        # Update main terrain visualization
        self.main_img.set_array(self.height_scaled[:, :, frame_idx])
        
        # Update title to show current generation step
        if frame_idx < len(self.descriptions):
            self.main_title.set_text(self.descriptions[frame_idx])
            
        # Parameters to visualize (raw and transformed)
        param_data_raw = [
            (self.raw_continentalness_scaled, self.continental_cmap),
            (self.raw_ridge_scaled, self.ridge_cmap),
            (self.raw_erosion_scaled, self.erosion_cmap),
            (self.raw_temperature_scaled, self.temp_cmap),
            (self.raw_humidity_scaled, self.humidity_cmap)
        ]
        
        param_data_transformed = [
            (self.continentalness_scaled, self.continental_cmap),
            (self.ridge_scaled, self.ridge_cmap),
            (self.erosion_scaled, self.erosion_cmap),
            (self.temperature_scaled, self.temp_cmap),
            (self.humidity_scaled, self.humidity_cmap)
        ]
        
        # Clear previous parameter images if they exist
        if not self.param_imgs_raw:
            # First time setup
            for i, (data, cmap) in enumerate(param_data_raw):
                ax = self.param_axes[i]
                # Split the axis into two subplots
                divider = ax.figure.add_gridspec(1, 2, wspace=0.05, hspace=0.05, left=ax.get_position().x0, 
                                                bottom=ax.get_position().y0, right=ax.get_position().x1, 
                                                top=ax.get_position().y1)
                
                ax_raw = self.fig.add_subplot(divider[0, 0])
                ax_raw.set_title("Raw", fontsize=8)
                ax_raw.axis('off')
                
                ax_transformed = self.fig.add_subplot(divider[0, 1])
                ax_transformed.set_title("Transformed", fontsize=8)
                ax_transformed.axis('off')
                
                # Create the images
                raw_img = ax_raw.imshow(data[:, :, frame_idx], cmap=cmap, vmin=0, vmax=1)
                transformed_img = ax_transformed.imshow(param_data_transformed[i][0][:, :, frame_idx], 
                                                       cmap=cmap, vmin=0, vmax=1)
                
                self.param_imgs_raw.append(raw_img)
                self.param_imgs_transformed.append(transformed_img)
        else:
            # Update existing images
            for i in range(len(self.param_imgs_raw)):
                raw_data = param_data_raw[i][0]
                transformed_data = param_data_transformed[i][0]
                
                self.param_imgs_raw[i].set_array(raw_data[:, :, frame_idx])
                self.param_imgs_transformed[i].set_array(transformed_data[:, :, frame_idx])
        
        # Redraw the figure
        self.fig.canvas.draw_idle()
        
    def next_frame(self, event):
        if self.current_frame < self.height_data.shape[2] - 1:
            self.update_frame(self.current_frame + 1)
            
    def prev_frame(self, event):
        if self.current_frame > 0:
            self.update_frame(self.current_frame - 1)
            
    def show(self):
        plt.show()

# Main function
if __name__ == "__main__":
    print("Creating transfer function terrain visualization...")
    viz = TransferFunctionVisualizer()
    viz.show()
    print("Visualization complete!")