# Load required packages
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(survival)
library(dplyr)
library(shinyjs)
library(readxl)
library(caret)
library(timeROC)

# ============================================================================
# Read data and build Cox model directly in Shiny
# ============================================================================

# Function to load and preprocess data
load_and_preprocess_data <- function() {
  # Read Excel data
  data <- read_excel("C:\\Users\\23899\\Desktop\\data_clean1(1).xlsx", sheet = 1)
  
  # Filter data: exclude samples with unknown outcome (value=2)
  all_data <- data %>%
    filter(`是否复发（0否1是2未知）` != 2,
           `是否死亡（0否1是2未知）` != 2)
  
  # Convert outcome and time variables to numeric
  all_data$`是否死亡（0否1是2未知）` <- as.numeric(all_data$`是否死亡（0否1是2未知）`)
  all_data$`OS(月)` <- as.numeric(all_data$`OS(月)`)
  
  # Define 9 high SHAP value features
  shap_features <- c(
    "放疗分组（0无1有2不详）",
    "年龄",
    "腋窝手术方式（0仅前哨1腋清2前+腋清）",
    "乳腺癌手术方式（1单切2保乳3改良根治4NSM）",
    "N分期",
    "患侧（1：左侧；2：右侧；3：双侧）",
    "病理分期(1:IA/IB;2:IIA/IIB/ⅢA/ⅢB/ⅢC;3:Ⅳ)",
    "家族史（无0，有1，2未知）",
    "Ki67数值（%）"
  )
  
  # Add survival time and outcome variables
  analysis_vars <- c(shap_features, "OS(月)", "是否死亡（0否1是2未知）")
  
  # Ensure data contains these variables
  all_data <- all_data[, intersect(colnames(all_data), analysis_vars)]
  
  # Variable type definitions
  categorical_vars <- c(
    "患侧（1：左侧；2：右侧；3：双侧）",
    "乳腺癌手术方式（1单切2保乳3改良根治4NSM）",
    "腋窝手术方式（0仅前哨1腋清2前+腋清）",
    "家族史（无0，有1，2未知）",
    "放疗分组（0无1有2不详）"
  )
  
  continuous_vars <- c(
    "年龄",
    "Ki67数值（%）",
    "N分期",
    "病理分期(1:IA/IB;2:IIA/IIB/ⅢA/ⅢB/ⅢC;3:Ⅳ)"
  )
  
  # Convert variable types
  for (var in categorical_vars) {
    if (var %in% colnames(all_data)) {
      all_data[[var]] <- factor(all_data[[var]])
    }
  }
  
  for (var in continuous_vars) {
    if (var %in% colnames(all_data)) {
      all_data[[var]] <- as.numeric(as.character(all_data[[var]]))
    }
  }
  
  # Remove missing values
  all_data <- na.omit(all_data)
  
  return(list(
    data = all_data,
    shap_features = shap_features,
    categorical_vars = categorical_vars,
    continuous_vars = continuous_vars
  ))
}

# ============================================================================
# Split data into training and test sets with 7:3 ratio (seed 43)
# ============================================================================

split_data <- function(all_data) {
  set.seed(43)
  
  # Split data with 7:3 ratio
  train_indices <- createDataPartition(all_data$`是否死亡（0否1是2未知）`, 
                                       p = 0.7, 
                                       list = FALSE)
  
  train_data <- all_data[train_indices, ]
  test_data <- all_data[-train_indices, ]
  
  return(list(train = train_data, test = test_data))
}

# ============================================================================
# Fit Cox model (using training set)
# ============================================================================

fit_cox_model <- function(train_data, shap_features) {
  
  # Build model formula
  processed_vars <- sapply(shap_features, function(x) {
    if (grepl("[^a-zA-Z0-9._]", x)) {
      return(paste0("`", x, "`"))
    } else {
      return(x)
    }
  })
  
  formula_shap <- paste("Surv(`OS(月)`, `是否死亡（0否1是2未知）`) ~", 
                        paste(processed_vars, collapse = " + "))
  
  # Fit Cox regression model
  cox_shap <- coxph(as.formula(formula_shap), data = train_data)
  
  return(cox_shap)
}

# ============================================================================
# Calculate baseline survival function (Breslow estimator) - using training set
# ============================================================================

calculate_baseline_survival <- function(cox_model, train_data) {
  
  # Create baseline patient data frame (all variables at minimum values)
  # Get factor levels
  breast_levels <- levels(train_data$`乳腺癌手术方式（1单切2保乳3改良根治4NSM）`)
  ax_levels <- levels(train_data$`腋窝手术方式（0仅前哨1腋清2前+腋清）`)
  side_levels <- levels(train_data$`患侧（1：左侧；2：右侧；3：双侧）`)
  family_levels <- levels(train_data$`家族史（无0，有1，2未知）`)
  radio_levels <- levels(train_data$`放疗分组（0无1有2不详）`)
  
  baseline_patient <- data.frame(
    `放疗分组（0无1有2不详）` = factor(0, levels = radio_levels),
    `年龄` = 1,
    `腋窝手术方式（0仅前哨1腋清2前+腋清）` = factor(0, levels = ax_levels),
    `乳腺癌手术方式（1单切2保乳3改良根治4NSM）` = factor(1, levels = breast_levels),
    `N分期` = 0,
    `患侧（1：左侧；2：右侧；3：双侧）` = factor(1, levels = side_levels),
    `病理分期(1:IA/IB;2:IIA/IIB/ⅢA/ⅢB/ⅢC;3:Ⅳ)` = 1,
    `家族史（无0，有1，2未知）` = factor(0, levels = family_levels),
    `Ki67数值（%）` = 1,
    check.names = FALSE
  )
  
  # Get baseline survival curve (Breslow estimator)
  baseline_surv <- survfit(cox_model, newdata = baseline_patient)
  
  # Extract survival probabilities at all event times
  baseline_df <- data.frame(
    time = baseline_surv$time,
    surv = baseline_surv$surv
  )
  
  # Add time 0
  baseline_df <- rbind(
    data.frame(time = 0, surv = 1),
    baseline_df
  )
  
  return(baseline_df)
}

# ============================================================================
# Calculate model performance metrics from training set (retained but not shown in UI)
# ============================================================================

calculate_model_performance <- function(cox_model, train_data, test_data) {
  
  # Training set C-index
  train_c_index <- summary(cox_model)$concordance[1]
  
  # Test set C-index
  test_concordance <- concordance(cox_model, newdata = test_data)
  test_c_index <- test_concordance$concordance
  
  # Calculate 3-year and 5-year AUC
  train_lp <- predict(cox_model, newdata = train_data, type = "lp")
  test_lp <- predict(cox_model, newdata = test_data, type = "lp")
  
  # Initialize AUC values
  auc_3yr_train <- NA
  auc_3yr_test <- NA
  auc_5yr_train <- NA
  auc_5yr_test <- NA
  
  # 3-year AUC
  tryCatch({
    train_roc_3 <- timeROC(
      T = train_data$`OS(月)`,
      delta = train_data$`是否死亡（0否1是2未知）`,
      marker = train_lp,
      cause = 1,
      times = 36,
      weighting = "marginal"
    )
    auc_3yr_train <- train_roc_3$AUC[2]
  }, error = function(e) {
    auc_3yr_train <- NA
  })
  
  tryCatch({
    test_roc_3 <- timeROC(
      T = test_data$`OS(月)`,
      delta = test_data$`是否死亡（0否1是2未知）`,
      marker = test_lp,
      cause = 1,
      times = 36,
      weighting = "marginal"
    )
    auc_3yr_test <- test_roc_3$AUC[2]
  }, error = function(e) {
    auc_3yr_test <- NA
  })
  
  # 5-year AUC
  tryCatch({
    train_roc_5 <- timeROC(
      T = train_data$`OS(月)`,
      delta = train_data$`是否死亡（0否1是2未知）`,
      marker = train_lp,
      cause = 1,
      times = 60,
      weighting = "marginal"
    )
    auc_5yr_train <- train_roc_5$AUC[2]
  }, error = function(e) {
    auc_5yr_train <- NA
  })
  
  tryCatch({
    test_roc_5 <- timeROC(
      T = test_data$`OS(月)`,
      delta = test_data$`是否死亡（0否1是2未知）`,
      marker = test_lp,
      cause = 1,
      times = 60,
      weighting = "marginal"
    )
    auc_5yr_test <- test_roc_5$AUC[2]
  }, error = function(e) {
    auc_5yr_test <- NA
  })
  
  return(list(
    train_c_index = train_c_index,
    test_c_index = test_c_index,
    auc_3yr_train = auc_3yr_train,
    auc_3yr_test = auc_3yr_test,
    auc_5yr_train = auc_5yr_train,
    auc_5yr_test = auc_5yr_test
  ))
}

# ============================================================================
# Extract coefficient table from Cox model
# ============================================================================

extract_coefficients <- function(cox_model) {
  summary_cox <- summary(cox_model)
  coef_table <- summary_cox$coefficients
  conf_int <- summary_cox$conf.int
  
  # Create coefficient data frame
  var_names <- rownames(coef_table)
  
  result <- data.frame(
    Variable = var_names,
    Coefficient = round(coef_table[, "coef"], 4),
    HR = round(conf_int[, "exp(coef)"], 4),
    CI_lower = round(conf_int[, "lower .95"], 4),
    CI_upper = round(conf_int[, "upper .95"], 4),
    p_value = round(coef_table[, "Pr(>|z|)"], 4)
  )
  
  return(result)
}

# ============================================================================
# Calculate individual survival probabilities
# ============================================================================

calculate_individual_survival <- function(hr, baseline_df) {
  
  # Cox model formula: S(t|X) = [S0(t)]^HR
  individual_df <- baseline_df
  individual_df$surv <- baseline_df$surv^hr
  individual_df$hr <- hr
  
  return(individual_df)
}

# ============================================================================
# Calculate hazard ratio based on input features (人工输入系数)
# ============================================================================

calculate_hr <- function(input) {
  risk_score <- 0
  
  # Age (relative to baseline 1)
  risk_score <- risk_score + (as.numeric(input$age) - 1) * 0.7595
  
  # N stage (relative to baseline 0)
  risk_score <- risk_score + (as.numeric(input$n_stage) - 0) * 0.4254
  
  # Pathological stage (relative to baseline 1)
  risk_score <- risk_score + (as.numeric(input$path_stage) - 1) * 0.9368
  
  # Ki67 (relative to baseline 1, 1 for <20%, 2 for >=20%)
  risk_score <- risk_score + (as.numeric(input$ki67) - 1) * 0.8398
  
  # Radiotherapy (relative to baseline 0)
  radio_val <- as.numeric(input$radiotherapy)
  if (radio_val == 1) {
    risk_score <- risk_score + 0.1625   # Yes
  } else if (radio_val == 2) {
    risk_score <- risk_score + 0.2372   # Unknown
  }
  # No (0) is baseline, coefficient = 0
  
  # Family history (relative to baseline 0)
  family_val <- as.numeric(input$family_history)
  if (family_val == 1) {
    risk_score <- risk_score + (-0.8402)   # Yes
  } else if (family_val == 2) {
    risk_score <- risk_score + 1.2290      # Unknown
  }
  # No (0) is baseline, coefficient = 0
  
  # Breast surgery - dummy variables (relative to baseline 1)
  br_val <- as.numeric(input$breast_surgery)
  if (br_val == 2) {
    risk_score <- risk_score + (-0.3494)   # Breast conserving
  } else if (br_val == 3) {
    risk_score <- risk_score + 0.0724      # Modified radical mastectomy (MRM)
  }
  # Simple mastectomy (1) is baseline, NSM (4) has no coefficient (assumed 0)
  
  # Axillary surgery - dummy variables (relative to baseline 0)
  ax_val <- as.numeric(input$ax_surgery)
  if (ax_val == 1) {
    risk_score <- risk_score + 0.2030      # ALND (Axillary Lymph Node Dissection)
  } else if (ax_val == 2) {
    risk_score <- risk_score + (-0.2498)   # Sentinel + ALND
  }
  # Sentinel only (0) is baseline, coefficient = 0
  
  # Side (relative to baseline 1)
  if (as.numeric(input$side) == 2) {
    risk_score <- risk_score + (-0.0356)   # Right side
  }
  # Left side (1) is baseline, coefficient = 0
  # Bilateral (3) has no coefficient (assumed 0)
  
  hr <- exp(risk_score)
  return(hr)
}

# ============================================================================
# Generate continuous monthly data for survival curve
# ============================================================================

generate_monthly_data <- function(step_data, max_time = 120) {
  
  # Generate all months
  all_times <- 0:max_time
  
  # For each month, find corresponding survival probability (step function rule)
  surv_at_times <- sapply(all_times, function(t) {
    valid_times <- step_data$time[step_data$time <= t]
    if(length(valid_times) == 0) return(1)
    last_time <- max(valid_times)
    step_data$surv[step_data$time == last_time]
  })
  
  data.frame(
    time = all_times,
    surv = surv_at_times
  )
}

# ============================================================================
# Variable coding definitions (for UI)
# ============================================================================

variable_coding <- list(
  age = list(
    name = "Age",
    levels = c("<40 years" = 1, 
               "40-59 years" = 2, 
               "60-74 years" = 3, 
               ">74 years" = 4),
    baseline = 1
  ),
  n_stage = list(
    name = "N stage",
    levels = c("N0" = 0, "N1" = 1, "N2" = 2, "N3" = 3),
    baseline = 0
  ),
  radiotherapy = list(
    name = "Radiotherapy",
    levels = c("No" = 0, "Yes" = 1, "Unknown" = 2),
    baseline = 0
  ),
  path_stage = list(
    name = "Pathological stage",
    levels = c("Stage I" = 1, "Stage II/III" = 2),
    baseline = 1
  ),
  ki67 = list(
    name = "Ki67 (%)",
    levels = c("Low (<20%)" = 1, "High (≥20%)" = 2),
    baseline = 1
  ),
  family_history = list(
    name = "Family history",
    levels = c("No" = 0, "Yes" = 1, "Unknown" = 2),
    baseline = 0
  ),
  breast_surgery = list(
    name = "Breast surgery",
    levels = c("Simple mastectomy" = 1, "Breast conserving" = 2, "Modified radical or NSM" = 3),
    baseline = 1
  ),
  ax_surgery = list(
    name = "Axillary surgery",
    levels = c("SLNB" = 0, "ALND" = 1, "SLNB+ALND" = 2),
    baseline = 0
  ),
  side = list(
    name = "Side",
    levels = c("Left" = 1, "Right" = 2),
    baseline = 1
  )
)

# ============================================================================
# UI Interface
# ============================================================================

ui <- fluidPage(
  theme = shinytheme("cerulean"),
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      .compact-select .selectize-control { margin-bottom: 5px; }
      .compact-select .form-group { margin-bottom: 8px; }
      .well { padding: 12px; margin-bottom: 10px; }
      h4 { margin-top: 5px; margin-bottom: 10px; }
    "))
  ),
  
  titlePanel("Triple-Negative Breast Cancer Cox Model Survival Prediction System"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      class = "compact-select",
      h4("Patient Characteristics", style = "color: #2c3e50; font-weight: bold; margin-top: 0;"),
      
      # First row: two inputs side by side
      fluidRow(
        column(6, selectInput("age", "Age Group", 
                              choices = variable_coding$age$levels,
                              selected = 1)),
        column(6, selectInput("n_stage", "N Stage",
                              choices = variable_coding$n_stage$levels,
                              selected = 0))
      ),
      
      # Second row: two inputs side by side
      fluidRow(
        column(6, selectInput("radiotherapy", "Radiotherapy",
                              choices = variable_coding$radiotherapy$levels,
                              selected = 0)),
        column(6, selectInput("path_stage", "Pathological Stage",
                              choices = variable_coding$path_stage$levels,
                              selected = 1))
      ),
      
      # Third row: two inputs side by side
      fluidRow(
        column(6, selectInput("ki67", "Ki67 (%)",
                              choices = variable_coding$ki67$levels,
                              selected = 1)),
        column(6, selectInput("family_history", "Family History",
                              choices = variable_coding$family_history$levels,
                              selected = 0))
      ),
      
      # Fourth row: two inputs side by side
      fluidRow(
        column(6, selectInput("breast_surgery", "Breast Surgery",
                              choices = variable_coding$breast_surgery$levels,
                              selected = 1)),
        column(6, selectInput("ax_surgery", "Axillary Surgery",
                              choices = variable_coding$ax_surgery$levels,
                              selected = 0))
      ),
      
      # Fifth row: Side alone
      fluidRow(
        column(6, selectInput("side", "Side",
                              choices = variable_coding$side$levels,
                              selected = 1))
      ),
      
      br(),
      actionButton("calculate", "Calculate Survival Probability", 
                   class = "btn-primary btn-lg btn-block",
                   style = "margin-top: 5px; margin-bottom: 0;")
    ),
    
    mainPanel(
      width = 8,
      
      # Survival curve tab only
      tabsetPanel(
        tabPanel("Survival Curve",
                 br(),
                 fluidRow(
                   column(12,
                          wellPanel(
                            plotOutput("surv_curve", height = "450px")
                          )
                   )
                 ),
                 fluidRow(
                   column(3,
                          wellPanel(
                            h4("1-Year Survival", style = "text-align: center; margin-top: 0; margin-bottom: 5px;"),
                            h3(textOutput("surv_1yr"), 
                               style = "color: #27ae60; font-weight: bold; text-align: center; font-size: 24px; margin-top: 0; margin-bottom: 0;")
                          )
                   ),
                   column(3,
                          wellPanel(
                            h4("3-Year Survival", style = "text-align: center; margin-top: 0; margin-bottom: 5px;"),
                            h3(textOutput("surv_3yr"), 
                               style = "color: #27ae60; font-weight: bold; text-align: center; font-size: 24px; margin-top: 0; margin-bottom: 0;")
                          )
                   ),
                   column(3,
                          wellPanel(
                            h4("5-Year Survival", style = "text-align: center; margin-top: 0; margin-bottom: 5px;"),
                            h3(textOutput("surv_5yr"), 
                               style = "color: #f39c12; font-weight: bold; text-align: center; font-size: 24px; margin-top: 0; margin-bottom: 0;")
                          )
                   ),
                   column(3,
                          wellPanel(
                            h4("10-Year Survival", style = "text-align: center; margin-top: 0; margin-bottom: 5px;"),
                            h3(textOutput("surv_10yr"), 
                               style = "color: #9b59b6; font-weight: bold; text-align: center; font-size: 24px; margin-top: 0; margin-bottom: 0;")
                          )
                   )
                 )
        )
      )
    )
  )
)

# ============================================================================
# Server logic
# ============================================================================

server <- function(input, output, session) {
  
  # Global reactive values, store model and data
  values <- reactiveValues(
    cox_model = NULL,
    baseline_survival = NULL,
    train_data = NULL,
    test_data = NULL,
    coef_table = NULL,
    performance = NULL,
    data_loaded = FALSE,
    loading_error = NULL
  )
  
  # Load data and fit model on startup
  observe({
    # Use tryCatch to capture errors
    tryCatch({
      # Load data
      data_list <- load_and_preprocess_data()
      all_data <- data_list$data
      shap_features <- data_list$shap_features
      
      # Split data into training and test sets (7:3, seed 43)
      split <- split_data(all_data)
      train_data <- split$train
      test_data <- split$test
      
      # Fit Cox model
      cox_model <- fit_cox_model(train_data, shap_features)
      
      # Calculate baseline survival function
      baseline_survival <- calculate_baseline_survival(cox_model, train_data)
      
      # Calculate model performance (retained but not shown in UI)
      performance <- calculate_model_performance(cox_model, train_data, test_data)
      
      # Extract coefficient table
      coef_table <- extract_coefficients(cox_model)
      
      # Store in reactive values
      values$cox_model <- cox_model
      values$baseline_survival <- baseline_survival
      values$train_data <- train_data
      values$test_data <- test_data
      values$coef_table <- coef_table
      values$performance <- performance
      values$data_loaded <- TRUE
      values$loading_error <- NULL
      
    }, error = function(e) {
      values$loading_error <- paste("Data loading error:", e$message)
      values$data_loaded <- FALSE
    })
  })
  
  # Display loading status
  output$loading_status <- renderText({
    if (is.null(values$loading_error)) {
      if (values$data_loaded) {
        return("Model loaded successfully!")
      } else {
        return("Loading model...")
      }
    } else {
      return(values$loading_error)
    }
  })
  
  # Calculate hazard ratio (使用人工输入系数)
  hr_val <- eventReactive(input$calculate, {
    calculate_hr(input)
  })
  
  # Calculate individual survival probabilities
  individual_surv <- reactive({
    req(values$baseline_survival, hr_val())
    calculate_individual_survival(hr_val(), values$baseline_survival)
  })
  
  # Generate monthly data
  individual_monthly <- reactive({
    req(individual_surv())
    generate_monthly_data(individual_surv(), max_time = 120)
  })
  
  baseline_monthly <- reactive({
    req(values$baseline_survival)
    generate_monthly_data(values$baseline_survival, max_time = 120)
  })
  
  # Output 1-year survival
  output$surv_1yr <- renderText({
    req(individual_monthly())
    surv <- individual_monthly()$surv[individual_monthly()$time == 12]
    if(length(surv) > 0) {
      paste0(round(surv * 100, 1), "%")
    } else {
      "N/A"
    }
  })
  
  # Output 3-year survival
  output$surv_3yr <- renderText({
    req(individual_monthly())
    surv <- individual_monthly()$surv[individual_monthly()$time == 36]
    if(length(surv) > 0) {
      paste0(round(surv * 100, 1), "%")
    } else {
      "N/A"
    }
  })
  
  # Output 5-year survival
  output$surv_5yr <- renderText({
    req(individual_monthly())
    surv <- individual_monthly()$surv[individual_monthly()$time == 60]
    if(length(surv) > 0) {
      paste0(round(surv * 100, 1), "%")
    } else {
      "N/A"
    }
  })
  
  # Output 10-year survival
  output$surv_10yr <- renderText({
    req(individual_monthly())
    surv <- individual_monthly()$surv[individual_monthly()$time == 120]
    if(length(surv) > 0) {
      paste0(round(surv * 100, 1), "%")
    } else {
      "N/A"
    }
  })
  
  # Plot survival curve
  output$surv_curve <- renderPlot({
    req(individual_monthly(), baseline_monthly(), hr_val(), values$baseline_survival)
    
    hr <- hr_val()
    indiv <- individual_monthly()
    baseline <- baseline_monthly()
    events <- values$baseline_survival
    
    # Key time points
    key_times <- c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120)
    
    # Get individual values at key time points
    indiv_key <- indiv[indiv$time %in% key_times, ]
    
    # Limit number of event points displayed to avoid overcrowding
    events_display <- events[seq(1, nrow(events), by = 3), ]
    
    ggplot() +
      # Individual survival curve
      geom_line(data = indiv, 
                aes(x = time, y = surv, color = "Individual patient"), 
                size = 1.2) +
      # Baseline survival curve
      geom_line(data = baseline, 
                aes(x = time, y = surv, color = "Baseline patient"), 
                size = 0.8, linetype = "dashed") +
      # Mark event points
      geom_point(data = events_display, 
                 aes(x = time, y = surv), 
                 color = "#7f8c8d", size = 1, alpha = 0.3) +
      # Mark key time points
      geom_point(data = indiv_key, 
                 aes(x = time, y = surv), 
                 color = "#2980b9", size = 3) +
      # Add value labels
      geom_text(data = indiv_key, 
                aes(x = time, y = surv, label = paste0(round(surv*100, 0), "%")), 
                vjust = -1, size = 3.5, color = "#2980b9") +
      scale_color_manual(values = c("Individual patient" = "#2980b9", "Baseline patient" = "#7f8c8d")) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                         limits = c(0, 1),
                         breaks = seq(0, 1, 0.1)) +
      scale_x_continuous(limits = c(0, 120), 
                         breaks = seq(0, 120, 12)) +
      labs(x = "Time (months)", y = "Survival Probability",
           title = "Individualized Survival Curve",
           subtitle = paste0("HR = ", round(hr, 3)),
           color = "") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray", size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
  })
  
  # Auto-calculate on page load
  observe({
    req(values$data_loaded)
    shinyjs::click("calculate")
  })
}

# ============================================================================
# Run the application
# ============================================================================

shinyApp(ui = ui, server = server)