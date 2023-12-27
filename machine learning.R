# Cargar bibliotecas
library(randomForest)
library(shiny)

# Cargar conjunto de datos iris
data(iris)

# Dividir el conjunto de datos en entrenamiento y prueba
set.seed(123)
training_indices <- sample(1:nrow(iris), 0.8 * nrow(iris))
train_data <- iris[training_indices, ]
test_data <- iris[-training_indices, ]

# Modelo de Random Forest
model <- randomForest(Species ~ ., data = train_data, ntree = 500)

# Shiny App
ui <- fluidPage(
  titlePanel("Predicción de la especie de Iris con Random Forest"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("input_sepal_length", "Longitud del sépalo", min = 4, max = 8, value = 5, step = 0.1),
      sliderInput("input_sepal_width", "Ancho del sépalo", min = 2, max = 4.5, value = 3, step = 0.1),
      sliderInput("input_petal_length", "Longitud del pétalo", min = 1, max = 7, value = 4, step = 0.1),
      sliderInput("input_petal_width", "Ancho del pétalo", min = 0.1, max = 2.5, value = 1, step = 0.1),
      actionButton("predict_button", "Predecir")
    ),
    mainPanel(
      textOutput("prediction_text")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$predict_button, {
    new_data <- data.frame(
      Sepal.Length = input$input_sepal_length,
      Sepal.Width = input$input_sepal_width,
      Petal.Length = input$input_petal_length,
      Petal.Width = input$input_petal_width
    )
    
    prediction <- predict(model, newdata = new_data, type = "response")
    
    output$prediction_text <- renderText({
      paste("Especie predicha:", prediction)
    })
  })
}

shinyApp(ui, server)
