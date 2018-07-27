使用教學:
可以直接執行hw4.exe。它會自動讀取同個路徑下名為input.txt的檔案當作輸入。若要自行指定檔名，可使用cmd   

 ./hw4 “檔名”   

 的指令，以參數的形式把要使用的input檔名傳入。(預設input.txt是複製rabbit.txt的內容)
(老師提供的兩個input我分別命名為monkey.txt和rabbit.txt)
程式執行完會輸出一個名為"colorOutput.ppm"的ppm圖片檔，請使用infraview來檢視。
專案中有一個codeblocks的專案檔hw4.cbp，可以使用此專案檔配合code blocks來開起來編譯。
Hw4.cpp : 這次作業的程式本體
algebra3.cpp: 這次作業使用老師提供的algebra3函式庫
imageio.cpp : ppm檔案相關的讀寫操作。

並且，我有稍微更改了input的格式
在每個"單一"物件的所有triangle mesh的敘述出現之後之後加上一個字元"O"的line來代表以上的所有MESH代表一個單一物件，需要有自己的KD-TREE。這樣做的目的是為了避免bounding box差距太大的兩個物體共用同一個kd-tree，會讓kd-tree的加速效果退化非常多。
例如老師的範例場景中，動物模型再靠右側，但是地板平面幾乎橫跨畫面的x方向，如果兩者共用同一個KD-tree，會變成動物模型任何一個kd-node中只要有包含到地板平面的mesh，他的bounding box幾乎就會橫跨整個畫面的x軸，造成很容易打中某些Bounding box，但是實際上該bounding box的mesh只有兩個三角形。

我把老師的monkey和bunny模型的模型本體跟地板的plane分別弄成兩個不同的物件
要測試請直接使用資料夾中的monkey.txt和rabbit.txt，裡面為我已經加上額外兩行"O"的物件描述的版本