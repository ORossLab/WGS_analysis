import argparse
import pandas as pd
import json
import matplotlib.pyplot as plt


# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-f", "--SV_file", help = "Path to SV file")

# Read arguments from command line
args = parser.parse_args()

print(args.SV_file)

chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
        'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

def remove_subset_groups(groups_dict):
    keys_to_remove = set()
    keys = list(groups_dict.keys())

    # Compare each group with every other group
    for i in range(len(keys)):
        for j in range(len(keys)):
            if i != j:
                # If group i is a subset of group j, mark it for removal
                if groups_dict[keys[i]].issubset(groups_dict[keys[j]]):
                    keys_to_remove.add(keys[i])

    # Remove the subset groups from the original dictionary
    filtered_groups = {key: groups_dict[key] for key in groups_dict if key not in keys_to_remove}

    return filtered_groups

def get_head_SVs(df):
    # Discover nested SVs:
    # Sort the DataFrame by start position
    df_sorted = df.sort_values(by='start').reset_index(drop=True)

    # # List to hold SVs that start and end before others
    before_sv_list = []

    # # Step 1: Compare each SV with all others
    for i in range(len(df_sorted)):
        sv_i = df_sorted.iloc[i]
        for j in range(len(df_sorted)):
            if i != j:  # Do not compare the same SV
                sv_j = df_sorted.iloc[j]
                # Check if sv_i starts and ends before sv_j
                if sv_i['start'] < sv_j['start'] and sv_i['end'] > sv_j['end']:
                    before_sv_list.append((sv_i['id'], sv_j['id']))  # Store the IDs of SVs


    # Step 2: Group nested SVs into sets
    nested_groups = {}

    for sv1, sv2 in before_sv_list:
        if sv1 not in nested_groups:
            nested_groups[sv1] = set([sv1])
        nested_groups[sv1].add(sv2)


    # Apply the function to remove subset groups
    filtered_nested_groups = remove_subset_groups(nested_groups)

    # Output the filtered groups
    if bool(filtered_nested_groups):
        print(f"Filtered Nested Groups for {chr}. Take a closer look to the following:")
        print(filtered_nested_groups)

        return filtered_nested_groups
    else:
        return {}

def convert_sets_to_lists(dict_obj):
    return {key: list(value) if isinstance(value, set) else value for key, value in dict_obj.items()}



vacmap_df = pd.read_csv(args.SV_file, sep="#")
vacmap_df = vacmap_df[(vacmap_df["chrom"].isin(chrs)) &
                      (vacmap_df["svlen"]>=10000)]

grouped_dfs = {chrom: group for chrom, group in vacmap_df.groupby('chrom')}
chrs = [value for value in chrs if value in grouped_dfs.keys()]


# Define the number of subplots (one for each chromosome)
num_chroms = len(grouped_dfs.keys())

# Create list to store nested SV info
final_nested_groups = []

# Create a figure and axes for subplots
fig, axes = plt.subplots(num_chroms, 1, figsize=(10, 2 * num_chroms), sharex=False)

# Loop through chromosomes and create a plot for each
for i, chr in enumerate(chrs):
    df = grouped_dfs[chr]
    filtered_nested_groups = get_head_SVs(df)
    final_nested_groups.append(filtered_nested_groups)

    # Create visualizations
    svtype_colors = {
        'DEL': 'red',
        'INS': 'blue',
        'INV': 'green',
        'DUP': 'purple',
        # Add other SV types as needed
    }

    num_bins = 100  # Adjust the number of bins

    # Create a list of start positions for each svtype
    svtypes = df['svtype'].unique()
    hist_data = [df[df['svtype'] == svtype]['start'] / 1_000_000 for svtype in svtypes]
    colors = [svtype_colors[svtype] for svtype in svtypes]

    # Plot stacked histograms on the appropriate axis
    axes[i].hist(hist_data, bins=num_bins, color=colors, stacked=True, edgecolor='white', label=svtypes, alpha=0.8)

    # Customize axis
    axes[i].set_title(f'{chr}')
    axes[i].set_ylabel('SV Count')
    axes[i].legend(title='SV Type')

    # Disable scientific notation for X-axis and label as millions
    axes[i].get_xaxis().get_major_formatter().set_scientific(False)



# Customize the X-axis for the bottom plot only
axes[-1].set_xlabel('Position on Chromosome (start) (x$\mathregular{10^{6}}$)')

# Apply a clean style to the entire figure
plt.style.use('ggplot')

# Adjust spacing between plots to prevent overlap
plt.tight_layout()

# Save the combined figure with all plots stacked
svg_file_path = "VACmap_SV_histogram.svg"
plt.savefig(svg_file_path, format='svg', bbox_inches='tight')  # 'bbox_inches' ensures no clipping of the plot elements
plt.close()

# Writing the list of dictionaries to a txt file=
with open("VACmap_nested_SVs.json", 'w') as f:
    for dictionary in final_nested_groups:
        # Convert the dictionary to a JSON string and write it to the file
        json_compatible_dict = convert_sets_to_lists(dictionary)
        f.write(json.dumps(json_compatible_dict) + '\n')
