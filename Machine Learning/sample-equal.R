

your_dataset <- read.csv(file = './data/data.csv')
desired_proportions <- c(SS = 0.25, S = 0.25, VS = 0.25, D = 0.25)

create_sample_equal_dataset <- function(your_dataset, desired_proportions) {
    #Extract the original proportions of C, S, and VS
    original_proportions <- table(your_dataset$Annotation) / nrow(your_dataset)
    #cat("- Original Proportions :")
    #original_proportions

    # Define the desired proportions for the bootstrap dataset
    #cat("\n- Desired Proportions :\n")
    #desired_proportions


    # Calculate the number of samples to draw for each group
    num_samples_per_group <- round(desired_proportions * nrow(your_dataset))
    #cat("\n- Number samples per group (desired ):\n")
    #num_samples_per_group


    # Initialize an empty dataframe for the bootstrap dataset
    bootstrap_dataset <- data.frame()

    # Loop through each group and sample the data
    for (group in names(num_samples_per_group)) {
        group_data <- your_dataset[your_dataset$Annotation == group, , drop = FALSE]
        sampled_data <- group_data[sample(1:nrow(group_data), num_samples_per_group[group], replace = TRUE), , drop = FALSE]
        bootstrap_dataset <- rbind(bootstrap_dataset, sampled_data)
    }

    # Shuffle the rows of the bootstrap dataset
    bootstrap_dataset <- bootstrap_dataset[sample(1:nrow(bootstrap_dataset)), ]

    new_proportions <- table(bootstrap_dataset$Annotation) / nrow(bootstrap_dataset)
    #cat("\n- New Proportions :")
    #new_proportions

    bootstrap_dataset
}

#create_sample_equal_dataset(your_dataset, desired_proportions)