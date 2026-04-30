window.onload = function () {
  loadTextFile("./calculate_star_temperature.txt", "calculate_code");
  loadTextFile("./download_data.txt", "download_code");
};

function loadTextFile(filePath, elementId) {
  fetch(filePath)
    .then(response => {
      if (!response.ok) {
        throw new Error("HTTP error: " + response.status);
      }
      return response.text();
    })
    .then(text => {
      document.getElementById(elementId).value = text;
    })
    .catch(error => {
      console.error(error);
      document.getElementById(elementId).value =
        filePath + " を読み込めませんでした。";
    });
}

function copyCode(elementId) {
  const targetCode = document.getElementById(elementId);

  navigator.clipboard.writeText(targetCode.value)
    .then(() => {
      alert("コピーしました！");
    })
    .catch(error => {
      console.error(error);
      alert("コピーに失敗しました。");
    });
}